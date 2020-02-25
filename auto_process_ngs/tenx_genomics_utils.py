#!/usr/bin/env python
#
#     tenx_genomics_utils.py: utility functions for handling 10xGenomics data
#     Copyright (C) University of Manchester 2017-2020 Peter Briggs
#

"""
tenx_genomics_utils.py

Utility classes and functions for processing the outputs from 10xGenomics's
Chromium SC 3'v2 system:

- MetricSummary
- AtacSummary
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
from builtins import range
from bcftbx.IlluminaData import SampleSheet
from bcftbx.IlluminaData import IlluminaData
from bcftbx.IlluminaData import IlluminaDataError
from bcftbx.IlluminaData import split_run_name_full
from bcftbx.JobRunner import SimpleJobRunner
from bcftbx.TabFile import TabFile
from bcftbx.utils import mkdirs
from bcftbx.utils import find_program
from .applications import Command
from .simple_scheduler import SchedulerJob
from .simple_scheduler import SimpleScheduler
from .simple_scheduler import SchedulerReporter
from .docwriter import Document
from .docwriter import List
from .docwriter import Link
from .docwriter import Table
from .analysis import AnalysisProject
from .bcl2fastq_utils import get_bases_mask
from .utils import get_numbered_subdir
from .utils import ZipArchive
from . import css_rules

# Initialise logging
import logging
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

# Permissible values for cellranger count --chemistry option
# See https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count
CELLRANGER_ASSAY_CONFIGS = {
    'auto': 'autodetection',
    'threeprime': 'Single Cell 3\'',
    'fiveprime': 'Single Cell 5\'',
    'SC3Pv1': 'Single Cell 3\' v1',
    'SC3Pv2': 'Single Cell 3\' v2',
    'SC3Pv3': 'Single Cell 3\' v3',
    'SC5P-PE': 'Single Cell 5\' paired-end (both R1 and R2 are used for alignment)',
    'SC5P-R2': 'Single Cell 5\' R2-only (where only R2 is used for alignment)',
}

#######################################################################
# Classes
#######################################################################

class MetricsSummary(object):
    """
    Extract data from metrics_summary.csv file

    Utility class for extracting data from a
    'metrics_summary.csv' file output from running
    'cellranger count'.

    The file consists of two lines: the first is a
    header line, the second consists of corresponding
    data values.

    In addition: integer data values are formatted to
    use commas to separate thousands (e.g. 2,272) and
    values which contain commas are enclosed in
    double quotes.

    For example:

    Estimated Number of Cells,Mean Reads per Cell,...
    "2,272","107,875","1,282","245,093,084",98.3%,...

    This class extracts the data values and where
    possible converts them to integers.
    """
    def __init__(self,s):
        """
        Create a new MetricsSummary instance

        Arguments:
          s (str): contents of a
            'metrics_summary.csv' file
        """
        self._data = dict()
        s = s.split('\n')
        fields = self._tokenise(s[0])
        line = self._tokenise(s[1])
        for key,value in zip(fields,line):
            self._data[key] = value
    def _tokenise(self,line):
        """
        Internal: process line from metrics_summary.csv

        Arguments:
          line (str): line to tokenise

        Returns:
          List: list of tokens extracted from the
            supplied line, with enclosing quotes
            removed and converted to integers
            where possible.
        """
        # Split the line into tokens
        # (taking account of quotes)
        tokens = list()
        this_token = ''
        quoted = False
        for c in line:
            if c == '"':
                # Start or end of a quote
                quoted = (not quoted)
            elif c == ',':
                if not quoted:
                    # End of a token
                    tokens.append(this_token)
                    this_token = ''
                    continue
            # Append to current token
            this_token += c
        # Deal with the last token
        tokens.append(this_token)
        # Strip quotes
        tokens = [t.strip('"') for t in tokens]
        # Convert to integer where possible
        # (i.e. remove commas from e.g. "2,272")
        for i in range(len(tokens)):
            try:
                tokens[i] = int(tokens[i].replace(',',''))
            except ValueError:
                pass
        return tokens
    @property
    def estimated_number_of_cells(self):
        """
        Return the estimated number of cells
        """
        return self._data['Estimated Number of Cells']

class AtacSummary(TabFile):
    """
    Extract data from summary.csv file for scATAC-seq

    Utility class for extracting data from a 'summary.csv'
    file output from running 'cellranger-atac count'.

    The file consists of two lines: the first is a
    header line, the second consists of corresponding
    data values.
    """
    def __init__(self,f):
        """
        Create a new AtacSummary instance

        Arguments:
          f (str): path to the 'summary.csv' file
        """
        TabFile.__init__(self,
                         filen=f,
                         first_line_is_header=True,
                         delimiter=',')
    @property
    def cells_detected(self):
        """
        Return the number of cells detected
        """
        return self[0]['cells_detected']
    @property
    def annotated_cells(self):
        """
        Return the number of annotated cells
        """
        return self[0]['annotated_cells']

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

def has_chromium_sc_indices(sample_sheet):
    """
    Check if a sample sheet contains Chromium SC indices

    The Chromium SC indices can be obtained from:

    https://support.10xgenomics.com/permalink/27rGqWvNYYuqkgeS66sksm

    The Chromium SC 3'v2 indices are of the form:

    SI-GA-[A-H][1-12]

    e.g. 'SI-GA-B11'

    For scATAC-seq the indices are of the form:

    SI-NA-[A-H][1-12]

    e.g. 'SI-NA-G9'

    Arguments:
      sample_sheet (str): path to the sample sheet CSV
        file to check

    Returns:
      Boolean: True if the sample sheet contains at least
        one Chromium SC index, False if not.
    """
    index_pattern = re.compile(r"SI-(GA|NA)-[A-H](1[0-2]|[1-9])$")
    s = SampleSheet(sample_sheet)
    for line in s:
        try:
            if index_pattern.match(line['index']):
                return True
        except KeyError:
            pass
    return False

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
        print(fastq)
        link_name = os.path.join(fq_dir,os.path.basename(fastq))
        if os.path.exists(link_name):
            logger.warning("%s: already exists" % link_name)
            continue
        target = os.path.relpath(fastq,fq_dir)
        logger.debug("Linking: %s -> %s" % (link_name,target))
        os.symlink(target,link_name)
    return fastq_path_dir

def set_cell_count_for_project(project_dir,qc_dir=None):
    """
    Set the total number of cells for a project

    Sums the number of cells for each sample in a project
    (as determined from 'cellranger count' and extracted
    from the 'metrics_summary.csv' file for scRNA-seq, or
    from 'cellranger-atac count' and extracted from the
    'summary.csv' file for scATAC-seq).

    The final count is written to the 'number_of_cells'
    metadata item for the project.

    Arguments:
      project_dir (str): path to the project directory
      qc_dir (str): path to QC directory (if not the default
        QC directory for the project)

    Returns:
      Integer: exit code, non-zero values indicate problems
        were encountered.
    """
    project = AnalysisProject(os.path.basename(project_dir),
                              project_dir)
    if qc_dir is None:
        qc_dir = project.qc_dir
    qc_dir = os.path.abspath(qc_dir)
    number_of_cells = 0
    if project.info.library_type in ('scRNA-seq',
                                     'snRNA-seq'):
        # Single cell/single nuclei RNA-seq
        for sample in project.samples:
            try:
                metrics_summary_csv = os.path.join(
                    qc_dir,
                    "cellranger_count",
                    sample.name,
                    "outs",
                    "metrics_summary.csv")
                metrics = MetricsSummary(
                    open(metrics_summary_csv,'r').read())
                number_of_cells += metrics.estimated_number_of_cells
            except Exception as ex:
                logger.critical("Failed to add cell count for sample "
                                "'%s': %s" % (sample.name,ex))
                return 1
    elif project.info.library_type in ('scATAC-seq',
                                       'snATAC-seq'):
        # Single cell ATAC-seq
        for sample in project.samples:
            try:
                summary_csv = os.path.join(
                    qc_dir,
                    "cellranger_count",
                    sample.name,
                    "outs",
                    "summary.csv")
                number_of_cells += AtacSummary(summary_csv).annotated_cells
            except Exception as ex:
                logger.critical("Failed to add cell count for sample "
                                "'%s': %s" % (sample.name,ex))
                return 1
    else:
        raise Exception("%s: don't know how to set cell count for "
                        "library type '%s'" % (project.name,
                                               project.info.library_type))
    # Store in the project metadata
    project.info['number_of_cells'] = number_of_cells
    project.info.save()
    return 0

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
            if line.startswith(package_name):
                # Extract version from line of the form
                # cellranger  (2.0.1)
                try:
                    package_version = line.split('(')[-1].strip(')')
                except Exception as ex:
                    logger.warning("Unable to get version from '%s': %s" %
                                   (line,ex))
    else:
        # No package supplied or located
        logger.warning("Unable to identify cellranger package "
                       "from '%s'" % cellranger_path)
    # Return what we found
    return (cellranger_path,package_name,package_version)

def run_cellranger_mkfastq(sample_sheet,
                           primary_data_dir,
                           output_dir,
                           lanes=None,
                           bases_mask=None,
                           ignore_dual_index=False,
                           minimum_trimmed_read_length=None,
                           mask_short_adapter_reads=None,
                           cellranger_exe='cellranger',
                           cellranger_jobmode='local',
                           cellranger_maxjobs=None,
                           cellranger_mempercore=None,
                           cellranger_jobinterval=None,
                           cellranger_localcores=None,
                           cellranger_localmem=None,
                           runner=None,
                           working_dir=None,
                           log_dir=None,
                           dry_run=False):
    """
    Wrapper for running 'cellranger mkfastq'

    Runs the 10xGenomics 'cellranger mkfastq' command to
    generate Fastqs from bcl files for Chromium single-cell
    data.

    To run the 'mkfastq' command using a different version
    of cellranger (e.g. cellranger-atac), specify the
    cellranger executable using the 'cellranger_exe'
    argument.

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
      ignore_dual_index (bool): optional, on a dual-indexed
        flowcell where the second index was not used for
        the 10x sample, ignore it
      minimum_trimmed_read_length: optional, specify minimum
        length for reads after adapter trimming (shorter reads
        will be padded with Ns to make them long enough)
      mask_short_adapter_reads: optional, specify the minimum
        length of ACGT bases that must be present in a read
        after adapter trimming for it not to be masked
        completely with Ns
      cellranger_exe (str): optional, name or path to
        cellranger executable (default: "cellranger")
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: "local")
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      cellranger_localcores (int): maximum number of cores
        cellranger can request in jobmode 'local'
        (default: None)
      cellranger_localmem (int): maximum memory cellranger
        can request in jobmode 'local' (default: None)
      runner (JobRunner): specify a non-default job runner
        to use when running cellranger (default:
        SimpleJobRunner)
      working_dir (str): path to a directory to use as
        as the working directory (default: current
        working directory)
      log_dir (str): path to a directory to write logs
        (default: current working directory)
      dry_run (bool): if True then only report actions
        that would be performed but don't run anything

    Returns:
      Integer: exit code from the cellranger command.
    """
    # Working directory
    if working_dir is None:
        working_dir = os.getcwd()
    # Job runner
    if runner is None:
        runner = SimpleJobRunner(join_logs=True)
    # Check for existing cellranger outputs
    flow_cell_dir = os.path.join(working_dir,
                                 flow_cell_id(primary_data_dir))
    if lanes is not None:
        lanes_suffix = "_%s" % lanes.replace(',','')
    else:
        lanes_suffix = ""
    flow_cell_dir = "%s%s" % (flow_cell_dir,lanes_suffix)
    mro_file = os.path.join(working_dir,
                            "__%s.mro" %
                            os.path.basename(flow_cell_dir))
    if not dry_run:
        if os.path.exists(flow_cell_dir):
            logger.warning("Removing existing output directory: %s" %
                           flow_cell_dir)
            shutil.rmtree(flow_cell_dir)
        if os.path.exists(mro_file):
            logger.warning("Removing existing mro file: %s" % mro_file)
            os.remove(mro_file)
    # Construct the cellranger command
    cmd = Command(cellranger_exe,
                  "mkfastq",
                  "--samplesheet",sample_sheet,
                  "--run",primary_data_dir,
                  "--output-dir",output_dir,
                  "--qc")
    if lanes is not None:
        cmd.add_args("--lanes=%s" % lanes)
    if bases_mask is not None:
        cmd.add_args("--use-bases-mask=%s" % bases_mask)
    if ignore_dual_index:
        cmd.add_args("--ignore-dual-index")
    if minimum_trimmed_read_length:
        cmd.add_args('--minimum-trimmed-read-length',
                     minimum_trimmed_read_length)
    if mask_short_adapter_reads:
        cmd.add_args('--mask-short-adapter-reads',
                     mask_short_adapter_reads)
    add_cellranger_args(cmd,
                        jobmode=cellranger_jobmode,
                        mempercore=cellranger_mempercore,
                        maxjobs=cellranger_maxjobs,
                        jobinterval=cellranger_jobinterval,
                        localcores=cellranger_localcores,
                        localmem=cellranger_localmem)
    # Run the command
    print("Running %s" % cmd)
    if not dry_run:
        # Sort out the working directory
        if working_dir is None:
            working_dir = os.getcwd()
        else:
            working_dir = os.path.abspath(working_dir)
        # Make a log directory
        if log_dir is None:
            log_dir = os.getcwd()
        else:
            log_dir = os.path.abspath(log_dir)
        # Submit the job
        cellranger_mkfastq_job = SchedulerJob(
            runner,
            cmd.command_line,
            name='cellranger_mkfastq',
            working_dir=working_dir,
            log_dir=log_dir)
        cellranger_mkfastq_job.start()
        try:
            cellranger_mkfastq_job.wait()
        except KeyboardInterrupt as ex:
            logger.warning("Keyboard interrupt, terminating cellranger")
            cellranger_mkfastq_job.terminate()
            raise ex
        exit_code = cellranger_mkfastq_job.exit_code
        print("cellranger mkfastq completed: exit code %s" % exit_code)
        if exit_code != 0:
            logger.error("cellranger mkfastq exited with an error")
            return exit_code
        # Check outputs and QC summary report
        if not os.path.isdir(flow_cell_dir):
            logger.error("No output directory '%s'" % flow_cell_dir)
            return -1
        json_file = os.path.join(flow_cell_dir,
                                 "outs",
                                 "qc_summary.json")
        if not os.path.exists(json_file):
            logger.error("cellranger mkfastq failed to make "
                         "JSON QC summary file (%s not found)"
                         % json_file)
            return -1
        # Make HTML QC summary
        html_file = os.path.join(working_dir,
                                 "cellranger_qc_summary%s.html" %
                                 lanes_suffix)
        if os.path.exists(html_file):
            logger.warning("Removing existing HTML QC summary file: %s"
                           % html_file)
            os.remove(html_file)
        make_qc_summary_html(json_file,html_file)
        if not os.path.exists(html_file):
            logger.error("Failed to create HTML QC summary file "
                         "(%s not found)" % html_file)
            return -1
        return exit_code

def run_cellranger_count(fastq_dir,
                         reference_data_path,
                         chemistry='auto',
                         cellranger_exe='cellranger',
                         cellranger_jobmode='local',
                         cellranger_maxjobs=None,
                         cellranger_mempercore=None,
                         cellranger_jobinterval=None,
                         cellranger_localcores=None,
                         cellranger_localmem=None,
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
      reference_data_path (str): path to the cellranger
        compatible transcriptome reference data
        directory (for scRNA-seq) or ATAC reference
        genome data (for scATAC-seq)
      chemistry (str): assay configuration (set to
        'auto' to let cellranger determine this
        automatically; ignored if not scRNA-seq)
      cellranger_exe (str): optional, name or path to
        cellranger executable (default: "cellranger")
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: "local")
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      cellranger_localcores (int): maximum number of cores
        cellranger can request in jobmode 'local'
        (default: None)
      cellranger_localmem (int): maximum memory cellranger
        can request in jobmode 'local' (default: None)
      max_jobs (int): maxiumum number of concurrent
        count jobs to run; also used for maximum number
        of jobs each count pipeline can run at once
        (default: 4)
      log_dir (str): path to a directory to write logs
        (default: current working directory)
      dry_run (bool): if True then only report actions
        that would be performed but don't run anything
      summary_only (bool): if True then only collect
        the output 'web_summary.html' and
        'metrics_summary.csv' files, otherwise
        copy all outputs (warning: this can be very
        large)

    Returns:
      Integer: exit code from the cellranger command.
    """
    # Cellranger mode
    cellranger_mode = os.path.basename(cellranger_exe)
    print("Cellranger mode: %s" % cellranger_mode)
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
        logger.critical("Couldn't load data from '%s'" %
                         fastq_dir)
        return 1
    print("Samples: %s" % sample_names)
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
        log_dir = get_numbered_subdir("%s_count" % cellranger_mode,
                                      parent_dir=log_dir,
                                      full_path=True)
        print("Log directory: %s" % log_dir)
        mkdirs(log_dir)

    # Submit the cellranger count jobs
    jobs = []
    for project in projects:
        print("Project: %s" % project)
        for sample in sample_names[project]:
            print("Sample: %s" % sample)
            # Check if outputs already exist
            count_dir = os.path.abspath(
                os.path.join(project,
                             "cellranger_count",
                             sample,
                             "outs"))
            if os.path.isdir(count_dir):
                print("-- %s: outputs exist, nothing to do" % sample)
                continue
            else:
                print("-- %s: setting up %s count" % (sample,
                                                      cellranger_mode))
            # Set up job for this sample
            work_dir = os.path.abspath("tmp.%s_count.%s.%s" %
                                       (cellranger_mode,
                                        project,sample))
            mkdirs(work_dir)
            print("Working directory: %s" % work_dir)
            cmd = Command(cellranger_exe,
                          "count",
                          "--id",sample,
                          "--fastqs",os.path.abspath(fastq_dir),
                          "--sample",sample)
            if cellranger_mode == "cellranger":
                cmd.add_args("--transcriptome",reference_data_path,
                             "--chemistry",chemistry)
            elif cellranger_mode == "cellranger-atac":
                cmd.add_args("--reference",reference_data_path)
            add_cellranger_args(cmd,
                                jobmode=cellranger_jobmode,
                                mempercore=cellranger_mempercore,
                                maxjobs=cellranger_maxjobs,
                                jobinterval=cellranger_jobinterval,
                                localcores=cellranger_localcores,
                                localmem=cellranger_localmem)
            print("Running: %s" % cmd)
            if not dry_run:
                job = sched.submit(cmd,
                                   name="%s_count.%s.%s" %
                                   (cellranger_mode,
                                    project,
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
        logger.critical("One or more jobs finished with non-zero "
                         "exit code")
        return retval

    # Handle outputs
    for project in projects:
        print("Project: %s" % project)
        for sample in sample_names[project]:
            print("Sample: %s" % sample)
            # Destination for count output
            count_dir = os.path.abspath(
                os.path.join(project,
                             "cellranger_count",
                             sample))
            mkdirs(count_dir)
            # Copy the cellranger count outputs
            outs_dir = os.path.join("tmp.%s_count.%s.%s"
                                    % (cellranger_mode,
                                       project,sample),
                                    sample,
                                    "outs")
            if not summary_only:
                # Collect all outputs
                print("Copying contents of %s to %s" % (outs_dir,count_dir))
                shutil.copytree(outs_dir,count_dir)
            else:
                # Only collect the web and csv summaries
                if cellranger_mode == 'cellranger':
                    files = ("web_summary.html","metrics_summary.csv")
                elif cellranger_mode == 'cellranger-atac':
                    files = ("web_summary.html","summary.csv")
                count_dir = os.path.join(count_dir,"outs")
                mkdirs(count_dir)
                for f in files:
                    path = os.path.join(outs_dir,f)
                    if not os.path.exists(path):
                        logger.warning("%s: not found in %s" % (f,outs_dir))
                        retval = 1
                    else:
                        print("Copying %s from %s to %s" % (f,
                                                            outs_dir,
                                                            count_dir))
                        shutil.copy(path,count_dir)
                # Stop if there was an error
                if retval != 0:
                    logger.critical("Some cellranger count outputs are "
                                    "missing")
                    return retval

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
        print("Making report for project %s" % project)
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
            print("Adding web summary (%s) for %s" % (web_summary,
                                                      sample))
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
                                     reference_data_path,
                                     chemistry='auto',
                                     cellranger_exe='cellranger',
                                     cellranger_jobmode='local',
                                     cellranger_maxjobs=None,
                                     cellranger_mempercore=None,
                                     cellranger_jobinterval=None,
                                     cellranger_localcores=None,
                                     cellranger_localmem=None,
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
      reference_data_path (str): path to the cellranger
        compatible transcriptome reference data
        directory (for scRNA-seq) or ATAC reference
        genome data (for scATAC-seq)
      chemistry (str): assay configuration (set to
        'auto' to let cellranger determine this
        automatically)
      cellranger_jobmode (str): specify the job mode to
        pass to cellranger (default: "local")
      cellranger_maxjobs (int): specify the maximum
        number of jobs to pass to cellranger (default:
        None)
      cellranger_mempercore (int): specify the memory
        per core (in Gb) to pass to cellranger (default:
        None)
      cellranger_jobinterval (int): specify the interval
        between launching jobs (in ms) to pass to
        cellranger (default: None)
      cellranger_localcores (int): maximum number of cores
        cellranger can request in jobmode 'local'
        (default: None)
      cellranger_localmem (int): maximum memory cellranger
        can request in jobmode 'local' (default: None)
      max_jobs (int): maxiumum number of concurrent
        count jobs to run; also used for maximum number
        of jobs each count pipeline can run at once
        (default: 4)
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
    retval = run_cellranger_count(
        fastq_dir,
        reference_data_path,
        chemistry=chemistry,
        cellranger_exe=cellranger_exe,
        cellranger_jobmode=cellranger_jobmode,
        cellranger_maxjobs=cellranger_maxjobs,
        cellranger_mempercore=cellranger_mempercore,
        cellranger_jobinterval=cellranger_jobinterval,
        cellranger_localcores=cellranger_localcores,
        cellranger_localmem=cellranger_localmem,
        max_jobs=max_jobs,
        log_dir=log_dir,
        dry_run=dry_run,
        summary_only=summary_only)
    # Extract and store the cell count from the cellranger
    # metric file
    if retval == 0:
        retval = set_cell_count_for_project(project_dir)
    return retval

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

import sys
if __name__ == "__main__":
    json_file = sys.argv[1]
    html_file = "qc_summary.html"
    report_cellranger_mkfastqs_qc_summary(json_file,html_file)
