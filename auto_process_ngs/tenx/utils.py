#!/usr/bin/env python
#
#     tenx/utils.py: utility functions for handling 10xGenomics data
#     Copyright (C) University of Manchester 2023-2024 Peter Briggs
#

"""
Utility functions for processing the outputs from 10x Genomics
pipelines:

- flow_cell_id
- has_10x_indices
- has_chromium_sc_indices
- get_bases_mask_10x_atac
- get_bases_mask_10x_multiome
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
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import split_run_name_full
from bcftbx.utils import find_program
from ..bcl2fastq.utils import get_bases_mask
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

def get_bases_mask_10x_atac(runinfo_xml):
    """
    Acquire a bases mask for 10xGenomics scATAC-seq

    Generates an initial bases mask based on the run
    contents, and then updates this so that:

    1. Only the first 8 bases of the first index read
       are actually used, and
    2. The second index read is converted to a data
       read.

    For example: if the initial bases mask is
    'Y50,I16,I16,Y50' then the scATAC-seq bases mask
    will be 'Y50,I8nnnnnnnn,Y16,Y50'.

    Arguments:
      runinfo_xml (str): path to the RunInfo.xml for
        the sequencing run

    Returns:
      String: 10xGenomics scATAC-seq bases mask string
    """
    bases_mask = get_bases_mask(runinfo_xml).lower().split(',')
    # Check there are four reads defined
    if len(bases_mask) != 4:
        raise Exception("Bases mask '%s' should have 4 reads "
                        "defined (has %d)" % (bases_mask,
                                              len(bases_mask)))
    # First read
    r1_mask = bases_mask[0]
    # Update first index to restrict to 8 bases
    num_cycles = int(bases_mask[1][1:])
    if num_cycles < 8:
        raise Exception("Index read < 8 bases")
    i1_mask = "I8%s" % ('n'*(num_cycles-8),)
    # Update second index to second read
    r2_mask = bases_mask[2].replace('i','y')
    # Keep last read as is
    r3_mask = bases_mask[3]
    # Reassemble and return
    return ','.join((r1_mask,i1_mask,r2_mask,r3_mask,))

def get_bases_mask_10x_multiome(runinfo_xml,library):
    """
    Return bases mask for 10xGenomics single cell multiome

    Generates an initial bases mask based on the run
    contents, and then updates this based on the library
    type (either 'atac' or 'gex').

    For ATAC data: the template bases mask is
    "Y*,I8n*,Y24,Y*" (keeping all of read 1, first 8 bases
    of read 2, all 24 bases of read 3, and all of read 4).

    For example: if the initial bases mask is
    'Y50,I10,Y24,Y90' then the single cell multiome ATAC
    bases mask will be 'Y50,I8n2,Y24,Y90'.

    For GEX data: the template bases mask is
    "Y28n*,I10,I10n*,Y*" (keeping first 28 bases of read 1,
    all 10 bases of read 2, first 10 bases of read 3, and
    all of read 4).

    For example: if the initial bases mask is
    'Y50,I10,Y24,Y90' then the single cell multiome GEX
    bases mask will be 'Y28n22,I10,I10n14,Y90'.

    Arguments:
      runinfo_xml (str): path to the RunInfo.xml for
        the sequencing run
      library (str): library type to set bases mask for
        (either 'atac' or 'gex')

    Returns:
      String: 10xGenomics single cell multiome bases mask
        string for the specified library type.
    """
    # Get initial bases mask from RunInfo.xml
    bases_mask = get_bases_mask(runinfo_xml).split(',')
    # Check there are four reads defined
    if len(bases_mask) != 4:
        raise Exception("Bases mask '%s' should have 4 reads "
                        "defined (has %d)" % (bases_mask,
                                              len(bases_mask)))
    # Regenerate bases mask
    library = library.lower()
    if library == "atac":
        # R1,R3: unchanged
        # R2: truncate to 24 bases
        # I1: truncate to 8 bases
        return get_bases_mask(runinfo_xml, i1=8, r2=24)
    elif library == "gex":
        # R1: truncate to 28 bases
        # I1,I2: truncate to 10 bases each
        # R2: unchanged
        return get_bases_mask(runinfo_xml, r1=28, i1=10, i2=10)
    else:
        raise Exception("Unknown library type: '%s'" % library)

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

def make_multi_config_template(f,reference=None,probe_set=None,
                               fastq_dir=None,samples=None,
                               no_bam=None,library_type="CellPlex",
                               cellranger_version=None):
    """
    Write a template configuration file for 'cellranger multi'

    Generates a template for the 'cellranger multi'
    configuration file, which can be used with either CellPlex
    or fixed RNA profiling (Flex) data.

    The format and parameters for different data types are
    described in the 10x Genomics 'cellranger' documentation:

    * CellPlex: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cellranger-multi
    * Flex: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi-frp#cellranger-multi

    Arguments:
      f (str): path that output template file will be
        written to
      reference (str): path to reference transcriptome
      probe_set (str): path to probe set CSV file
      fastq_dir (str): path to directory with Fastq files
      samples (list): list of sample names
      no_bam (bool): if set then will be the value of the
        'no-bam' setting
      library_type (str): specify the library type of data
        that the configuration file will be used with;
        should be one of 'CellPlex[...]' (the default),
        'Flex' or 'Single Cell Immune Profiling'
      cellranger_version (str): optionally specify the
        target Cellranger version number (or None)
    """
    # Target version
    if cellranger_version is None:
        cellranger_version = DEFAULT_CELLRANGER_VERSION
    else:
        cellranger_version = str(cellranger_version)
    # Major version
    try:
        cellranger_major_version = int(cellranger_version.split(".")[0])
    except Exception as ex:
        raise Exception(f"'{cellranger_version}': unable to extract "
                        f"Cellranger major  version from version string: "
                        f"{ex}")
    # Normalise and check supplied library type
    if library_type.startswith("CellPlex"):
        library_type = "cellplex"
    library_type = library_type.lower()
    if library_type not in ("cellplex",
                            "flex",
                            "single cell immune profiling"):
        raise Exception(f"'{library_type}': unsupported library type for "
                        "multi template generation")
    # Write template file
    with open(f,'wt') as fp:
        # Header
        fp.write("## 10x_multi_config.csv\n"
                 "## See:\n")
        if library_type == "cellplex":
                 fp.write("## * CellPlex: https://support.10xgenomics.com/"
                          "single-cell-gene-expression/software/pipelines/"
                          "latest/using/multi#cellranger-multi\n")
        elif library_type == "flex":
            fp.write("## * Flex: https://support.10xgenomics.com/"
                     "single-cell-gene-expression/software/pipelines/"
                     "latest/using/multi-frp#cellranger-multi\n")
        elif library_type == "single cell immune profiling":
            fp.write("## * Single Cell Immune Profiling: "
                     "https://support.10xgenomics.com/"
                     "single-cell-vdj/software/pipelines/latest/using/"
                     "multi\n")
        # Gene expression section
        fp.write("[gene-expression]\n")
        fp.write("reference,%s\n" %
                 (reference if reference
                  else "/path/to/transcriptome"))
        if library_type == "flex":
            fp.write("probe-set,%s\n" %
                     (probe_set if probe_set
                      else "/path/to/probe/set"))
        fp.write("#force-cells,n\n")
        if cellranger_major_version == 7:
            # Cellranger 7.* targetted
            if no_bam is not None:
                fp.write("no-bam,%s\n" % str(no_bam).lower())
            else:
                fp.write("#no-bam,true|false\n")
        elif cellranger_major_version in (8, 9):
            # Cellranger 8.* or 9.* targetted
            if no_bam is not None:
                fp.write("create-bam,%s\n" % str(not no_bam).lower())
            else:
                fp.write("create-bam,true\n")
        fp.write("#cmo-set,/path/to/custom/cmo/reference\n")
        fp.write("\n")
        # Feature section
        fp.write("#[feature]\n"
                 "#reference,/path/to/feature/reference\n")
        fp.write("\n")
        # V(D)J section
        if library_type == "single cell immune profiling":
            fp.write("#[vdj]\n"
                     "#reference,/path/to/vdj/reference\n")
            fp.write("\n")
        # Libraries section
        fp.write("[libraries]\n"
                 "fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate\n")
        if samples:
            if library_type == "cellplex":
                tenx_library_type = "[Gene Expression|Multiplexing Capture]"
            elif library_type == "flex":
                tenx_library_type = "[Gene Expression|Antibody Capture]"
            elif library_type == "single cell immune profiling":
                tenx_library_type = "[Gene Expression|Antibody Capture|VDJ-B|VDJ-T]"
            for sample in samples:
                fp.write("{sample},{fastqs_dir},any,{sample},{tenx_library_type},\n".format(
                    sample=sample,
                    fastqs_dir=(fastq_dir if fastq_dir else "/path/to/fastqs"),
                    tenx_library_type=tenx_library_type))
        # Multiplexed samples section
        if library_type == "cellplex":
            fp.write("\n[samples]\n")
            fp.write("sample_id,cmo_ids,description\n"
                     "MULTIPLEXED_SAMPLE,CMO1|CMO2|...,DESCRIPTION\n")
        elif library_type == "flex":
            fp.write("\n[samples]\n")
            fp.write("sample_id,probe_barcode_ids,description\n"
                     "MULTIPLEXED_SAMPLE,BC001|BC002|...,DESCRIPTION\n")
