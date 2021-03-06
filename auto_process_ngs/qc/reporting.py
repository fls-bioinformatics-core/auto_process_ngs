#!/usr/bin/env python
#
#     reporting: report QC from analysis projects
#     Copyright (C) University of Manchester 2018-2021 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import logging
import time
import ast
import shutil
from collections import defaultdict
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.TabFile import TabFile
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from ..analysis import AnalysisFastq
from ..analysis import run_reference_id
from ..analysis import split_sample_name
from ..docwriter import Document
from ..docwriter import Section
from ..docwriter import Table
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import Target
from ..docwriter import List
from ..docwriter import Para
from ..docwriter import WarningIcon
from ..metadata import AnalysisDirMetadata
from ..metadata import AnalysisProjectQCDirInfo
from ..fastq_utils import group_fastqs_by_name
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .fastq_strand import Fastqstrand
from .cellranger import CellrangerCount
from .outputs import fastqc_output
from .outputs import fastq_screen_output
from .outputs import fastq_strand_output
from .outputs import expected_outputs
from .plots import uscreenplot
from .plots import ufastqcplot
from .plots import uboxplot
from .plots import ustrandplot
from .plots import encode_png
from ..utils import ZipArchive
from .. import get_version

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Data
#######################################################################

from .constants import FASTQ_SCREENS

QC_REPORT_CSS_STYLES = """/* General styles */
html { font-family: DejaVu Serif, serif; }
p { font-size: 85%;
    color: #808080; }
/* Headers */
h1, h2, h3, h4 { font-family: DejaVu Serif, serif; }
h1 { background-color: #42AEC2;
     color: white;\n
     padding: 5px 10px; }
h2 { background-color: #8CC63F;
     color: white;
     display: inline-block;
     padding: 5px 15px;
     margin: 0;
     border-top-left-radius: 20px;
     border-bottom-right-radius: 20px; }
h3, h4 { background-color: grey;
         color: white;
         display: block;
         padding: 5px 15px;
         margin: 0;
         border-top-left-radius: 20px;
         border-bottom-right-radius: 20px; }
/* Summary section */
div.summary { margin: 10 10;
              border: solid 2px #8CC63F;
              padding: 0;
              border-top-left-radius: 25px; }
.info { padding: 5px 15px;
        float: left;
        font-size: 100%; }
.info h3 { margin: 5px; }
.info p { color: black; }
/* Samples and Fastqs */
.sample { margin: 10 10;
          border: solid 2px #8CC63F;
          padding: 0;
          background-color: #ffe;
          border-top-left-radius: 25px;
          border-bottom-right-radius: 25px; }
.fastqs { border: 1px solid grey;
          padding: 5px;
          margin: 5px 20px; }
.fastq { border: 2px solid lightgray;
         padding: 5px;
         margin: 5px;
         float: left; }
.strandedness { border: 2px solid lightgray;
                padding: 5px;
                margin: 5px;float: left; }
/* Metadata table */
table.metadata {
          margin: 10 10;
          border: solid 1px grey;
          background-color: white;
          font-size: 90%; }
table.metadata tr td:first-child {
          background-color: grey;
          color: white;
          padding: 2px 5px;
          font-weight: bold; }
/* Summary table */
table.summary { border: solid 1px grey;
          background-color: white;
          font-size: 80%; }
table.summary th { background-color: grey;
                   color: white;
                   padding: 2px 5px; }
table.summary td { text-align: center;
                   padding: 2px 5px;
                   border-bottom: solid 1px lightgray; }
table.summary tr td:first-child { text-align: right; }
table.fastq_summary tr td:first-child {
          background-color: grey;
          color: white;
          font-weight: bold; }
table.fastq_summary tr td:first-child a {
          color: white;
          font-weight: bold; }
/* Warnings section */
.warnings { padding: 2px;
            border: solid 3px red;
            color: red;
            background-color: #F5BCA9;
            font-weight: bold;
            margin: 10px; }
.warnings p   { color: red;
                font-size: 120%; }
.warnings img { vertical-align: middle; }
/* Display control elements */
.clear { clear: both; }
.hide  { display: none; }
/* FastQC summary table */
table.fastqc_summary span.PASS { font-weight: bold;
                                 color: green; }
table.fastqc_summary span.WARN { font-weight: bold;
                                 color: orange; }
table.fastqc_summary span.FAIL { font-weight: bold;
                                 color: red; }
/* Program versions */
table.programs th { text-align: left;
                    background-color: grey;
                    color: white;
                    padding: 2px 5px; }
table.programs td { padding: 2px 5px;
                    border-bottom: solid 1px lightgray; }
/* Rules for printing */
@media print
{
a { color: black; text-decoration: none; }
.sample { page-break-before: always; }
table th { border-bottom: solid 1px lightgray; }
.no_print { display: none; }
}
"""

#######################################################################
# Classes
#######################################################################

class QCReporter(object):
    """
    Class describing QC results for an AnalysisProject

    Provides the follow properties:

    name: project name
    paired_end: True if project is paired-end
    samples: list of sample names

    Provides the following methods:

    verify: checks the QC outputs for the project
    report: generate a HTML report for the project
    """
    def __init__(self,project):
        """
        Initialise a new QCReporter instance

        Arguments:
           project (AnalysisProject): project to report QC for
        """
        self._project = project

    def verify(self,qc_dir=None,qc_protocol=None):
        """
        Check the QC outputs are correct for the project

        Arguments:
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project being checked.
          qc_protocol (str): QC protocol to verify against
            (optional)

        Returns:
          Boolean: Returns True if all expected QC products
            are present, False if not.
        """
        return verify(self._project,
                      qc_dir=qc_dir,
                      qc_protocol=qc_protocol)

    def report(self,title=None,filename=None,qc_dir=None,
               report_attrs=None,summary_fields=None,
               relative_links=False,make_zip=False):
        """
        Report the QC for the project

        Arguments:
          title (str): optional, specify title for the report
            (defaults to '<PROJECT_NAME>: QC report')
          filename (str): optional, specify path and name for
            the output report file (defaults to
            '<PROJECT_NAME>.qc_report.html')
          qc_dir (str): path to the QC output dir
          report_attrs (list): optional, list of elements to
            report for each Fastq pair
          summary_fields (list): optional, list of fields to
            report for each sample in the summary table
          relative_links (boolean): optional, if set to True
            then use relative paths for links in the report
            (default is to use absolute paths)
          make_zip (boolean): if True then also create a ZIP
            archive of the QC report and outputs (default is
            not to create the ZIP archive)

        Returns:
          String: filename of the output HTML report.
        """
        return report((self._project,),
                      title=title,
                      filename=filename,
                      qc_dir=qc_dir,
                      report_attrs=report_attrs,
                      summary_fields=summary_fields,
                      relative_links=relative_links,
                      make_zip=make_zip)

class QCProject(object):
    """
    Gather information about the QC associated with a project

    Collects data about the QC for an AnalysisProject and
    makes it available via the following properties:

    - project: the AnalysisProject instance
    - qc_dir: the directory to examine for QC outputs
    - run_metadata: AnalysisDirMetadata instance with
      metadata from the parent run (if present)
    - processing_software: dictionary with information on
      software used to process the initial data

    Properties that shortcut to properties of the parent
    AnalysisProject:

    - name: project name
    - dirn: path to associated directory
    - comments: comments associated with the project
    - info: shortcut to the project's AnalysisProjectMetadata
      instance
    - qc_info: shortcut to the QC directory's
      AnalysisProjectQCDirInfo instance
    - fastq_attrs: class to use for extracting data from
      Fastq file names

    Properties based on artefacts detected in the QC
    directory:

    - fastqs: sorted list of Fastq names
    - reads: list of reads (e.g. 'r1', 'r2', 'i1' etc)
    - samples: sorted list of sample names extracted
      from Fastqs
    - outputs: list of QC output categories detected (see
      below for valid values)
    - output_files: list of absolute paths to QC output
      files
    - software: dictionary with information on the
      QC software packages

    The following are valid values for the 'outputs'
    property:

    - 'fastqc_[r1...]'
    - 'screens_[r1...]'
    - 'strandedness'
    - 'icell8_stats'
    - 'icell8_report'
    - 'cellranger_count'
    - 'multiqc'

    General properties about the project:

    - is_single_cell: True if the project has single cell
        data (10xGenomics, ICELL8 etc)
    """
    def __init__(self,project,qc_dir=None):
        """
        Create a new QCProject instance

        Arguments:
          project (AnalysisProject): project to report QC for
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
        """
        logger.debug("QCProject: project         : %s" % project.name)
        logger.debug("QCProject: project dir     : %s" % project.dirn)
        logger.debug("QCProject: qc_dir (initial): %s" % qc_dir)
        # Store project
        self.project = project
        # Sort out target QC dir
        if qc_dir is None:
            qc_dir = self.project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(self.project.dirn,
                                      qc_dir)
        self.qc_dir = qc_dir
        logger.debug("QCProject: qc_dir (final): %s" % self.qc_dir)
        # How to handle Fastq names
        self.fastq_attrs = self.project.fastq_attrs
        # Detect outputs
        self._detect_outputs()
        # Expose metadata from parent run
        self.run_metadata = AnalysisDirMetadata()
        run_metadata_file = os.path.join(
            os.path.dirname(os.path.abspath(self.project.dirn)),
            "metadata.info")
        if os.path.exists(run_metadata_file):
            print("Loading run metadata from %s" % run_metadata_file)
            self.run_metadata.load(run_metadata_file)
        else:
            logger.warning("Run metadata file '%s' not found"
                           % run_metadata_file)
        # Collect processing software metadata
        try:
            self.processing_software = ast.literal_eval(
                self.run_metadata.processing_software)
        except ValueError:
            self.processing_software = dict()
        if not self.processing_software:
            # Fallback to legacy metadata items
            try:
                self.processing_software['bcl2fastq'] = ast.literal_eval(
                    self.run_metadata.bcl2fastq_software)
            except ValueError:
                pass
            try:
                self.processing_software['cellranger'] = ast.literal_eval(
                    self.run_metadata.cellranger_software)
            except ValueError:
                pass
        # Expose project metadata
        self.info = self.project.info
        # Expose QC info
        self.qc_info = self.project.qc_info(self.qc_dir)

    @property
    def id(self):
        """
        Identifier for the project

        This is a string of the form:

        <RUN_ID>:<PROJECT_NAME>

        e.g. ``MINISEQ_201120#22:PJB``

        If the run id can't be determined then the name
        of the parent directory is used instead.
        """
        run_id = self.run_id
        if not run_id:
            run_id = os.path.basename(os.path.dirname(self.dirn))
        return "%s:%s" % (run_id,self.name)

    @property
    def name(self):
        """
        Name of project
        """
        return self.project.name

    @property
    def dirn(self):
        """
        Path to project directory
        """
        return self.project.dirn

    @property
    def comments(self):
        """
        Comments associated with the project
        """
        return self.project.info.comments

    @property
    def run_id(self):
        """
        Identifier for parent run

        This is the standard identifier constructed
        from the platform, datestamp and facility
        run number (e.g. ``MINISEQ_201120#22``).

        If an identifier can't be constructed then
        ``None`` is returned.
        """
        try:
            return run_reference_id(self.info['run'],
                                    platform=self.info['platform'],
                                    facility_run_number=
                                    self.run_metadata['run_number'])
        except (AttributeError,TypeError) as ex:
            logger.warning("Run reference ID can't be "
                           "determined: %s (ignored)" % ex)
            return None

    @property
    def is_single_cell(self):
        """
        Check whether project has single cell data
        """
        return (self.info['single_cell_platform'] is not None)

    def software_info(self,pkg):
        """
        Get information on software package

        Arguments:
          pkg (string): name of software package to
            get information about

        Returns:
          String: software version information, or
            None if no information is stored.
        """
        # Acquire the value
        try:
            if self.software[pkg]:
                return ','.join(self.software[pkg])
        except KeyError:
            try:
                return self.processing_software[pkg][2]
            except KeyError:
                # Not processing software
                return None
            except TypeError:
                # Not valid data
                return None

    def _detect_outputs(self):
        """
        Internal: determine which QC outputs are present
        """
        outputs = set()
        output_files = []
        software = {}
        print("Scanning contents of %s" % self.qc_dir)
        files = [os.path.join(self.qc_dir,f)
                 for f in os.listdir(self.qc_dir)]
        print("\t- %d objects found" % len(files))
        logger.debug("files: %s" % files)
        # Look for screen files
        screens = list(filter(lambda f:
                              f.endswith("_screen.txt") or
                              f.endswith("_screen.png"),
                              files))
        logger.debug("Screens: %s" % screens)
        print("\t- %d fastq_screen files" % len(screens))
        fastq_names = set()
        if screens:
            versions = set()
            # Pull out the Fastq names from the .txt files
            for screen in list(filter(lambda s:
                                      s.endswith("_screen.txt"),
                                      screens)):
                screen_base = os.path.splitext(screen)[0]
                s = os.path.basename(screen_base)[:-len("_screen")]
                for name in FASTQ_SCREENS:
                    if s.endswith("_%s" % name):
                        fq = self.fastq_attrs(s[:-len("_%s" % name)])
                        outputs.add("screens_%s%s" %
                                    (('i' if fq.is_index_read else 'r'),
                                     (fq.read_number
                                      if fq.read_number is not None else '1')))
                        fastq_names.add(s[:-len("_%s" % name)])
                versions.add(Fastqscreen(screen).version)
            if versions:
                software['fastq_screen'] = sorted(list(versions))
            # Store the screen files
            output_files.extend(screens)
        # Look for fastqc outputs
        fastqcs = list(filter(lambda f:
                              f.endswith("_fastqc") and
                              os.path.exists("%s.html" % f) and
                              os.path.exists(os.path.join(f,"summary.txt")) and
                              os.path.exists(os.path.join(f,"fastqc_data.txt")),
                              files))
        logger.debug("Fastqc: %s" % fastqcs)
        print("\t- %d fastqc files" % len(fastqcs))
        if fastqcs:
            versions = set()
            # Pull out the Fastq names from the Fastqc files
            for fastqc in fastqcs:
                f = os.path.basename(fastqc)[:-len("_fastqc")]
                fastq_names.add(f)
                fq = self.fastq_attrs(f)
                outputs.add("fastqc_%s%s" %
                            (('i' if fq.is_index_read else 'r'),
                             (fq.read_number
                              if fq.read_number is not None else '1')))
                versions.add(Fastqc(fastqc).version)
            if versions:
                software['fastqc'] = sorted(list(versions))
            # Store the fastqc files needed for reporting
            output_files.extend(["%s.html" % f for f in fastqcs])
            output_files.extend([os.path.join(f,"summary.txt")
                                 for f in fastqcs])
            output_files.extend([os.path.join(f,"fastqc_data.txt")
                                 for f in fastqcs])
            # Fastqc plot images
            for png in ("per_base_quality",):
                output_files.extend([os.path.join(f,
                                                  "Images",
                                                  "%s.png" % png)
                                     for f in fastqcs])
        # Look for fastq_strand outputs
        fastq_strand = list(filter(lambda f:
                                   f.endswith("_fastq_strand.txt"),
                                   files))
        logger.debug("fastq_strand: %s" % fastq_strand)
        print("\t- %d fastq_strand files" % len(fastq_strand))
        if fastq_strand:
            outputs.add("strandedness")
            versions = set()
            for f in fastq_strand:
                fq = self.fastq_attrs(os.path.splitext(f)[0])
                fastq_names.add(
                    os.path.basename(
                        os.path.splitext(f)[0])[:-len("_fastq_strand")])
                versions.add(Fastqstrand(f).version)
            if versions:
                software['fastq_strand'] = sorted(list(versions))
            # Store the fastq_strand files
            output_files.extend(fastq_strand)
        # Look for ICELL8 outputs
        print("Checking for ICELL8 reports in %s/stats" %
              self.project.dirn)
        icell8_stats_xlsx = os.path.join(self.project.dirn,
                                         "stats",
                                         "icell8_stats.xlsx")
        if os.path.exists(icell8_stats_xlsx):
            outputs.add("icell8_stats")
            output_files.append(icell8_stats_xlsx)
        icell8_report_html = os.path.join(self.project.dirn,
                                          "icell8_processing.html")
        if os.path.exists(icell8_report_html):
            outputs.add("icell8_report")
            output_files.append(icell8_report_html)
        # Look for cellranger_count outputs
        cellranger_count_dir = os.path.join(self.qc_dir,
                                            "cellranger_count")
        cellranger_samples = []
        if os.path.isdir(cellranger_count_dir):
            for d in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_count_dir,f)),
                    os.listdir(cellranger_count_dir)):
                web_summary_html = os.path.join(cellranger_count_dir,
                                                d,
                                                "outs",
                                                "web_summary.html")
                if os.path.exists(web_summary_html):
                    cellranger_samples.append(d)
                    output_files.append(web_summary_html)
            if cellranger_samples:
                outputs.add("cellranger_count")
        # Look for MultiQC report
        print("Checking for MultiQC report in %s" % self.project.dirn)
        multiqc_report = os.path.join(self.project.dirn,
                                      "multi%s_report.html"
                                      % os.path.basename(self.qc_dir))
        if os.path.isfile(multiqc_report):
            outputs.add("multiqc")
            output_files.append(multiqc_report)
        # Fastqs sorted by sample name
        self.fastqs = sorted(list(fastq_names),
                             key=lambda fq: split_sample_name(
                                 self.fastq_attrs(fq).sample_name))
        # Determine reads
        reads = set()
        for fastq in self.fastqs:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                reads.add("i%s" % (fq.read_number
                                   if fq.read_number is not None else '1'))
            else:
                reads.add("r%s" % (fq.read_number
                                   if fq.read_number is not None else '1'))
        self.reads = sorted(list(reads))
        # Samples
        samples = set([self.fastq_attrs(fq).sample_name
                       for fq in self.fastqs] +
                      [s.name for s in self.project.samples])
        for s in cellranger_samples:
            samples.add(s)
        self.samples = sorted(list(samples),
                              key=lambda s: split_sample_name(s))
        # QC outputs
        self.outputs = sorted(list(outputs))
        # Software versions
        self.software = software
        # Output files
        self.output_files = sorted(output_files)

class FastqSet(object):
    """
    Class describing a set of Fastq files

    A set can be a single or a pair of fastq files.

    Provides the following properties:

    r1: R1 Fastq in the pair
    r2: R2 Fastq (will be None if no R2)
    fastqs: list of Fastq files

    Provides the following methods:

    verify: checks the QC outputs for the set
    """
    def __init__(self,fqr1,fqr2=None):
        """
        Initialise a new QCFastqSet instance

        Arguments:
           fqr1 (str): path to R1 Fastq file
           fqr2 (str): path to R2 Fastq file, or
             None if the 'set' is a single Fastq
        """
        self._fastqs = list((fqr1,fqr2))

    def __getitem__(self,key):
        return self.fastqs[key]

    @property
    def r1(self):
        """
        Return R1 Fastq file from pair
        """
        return self._fastqs[0]

    @property
    def r2(self):
        """
        Return R2 Fastq file from pair
        """
        return self._fastqs[1]

    @property
    def fastqs(self):
        """
        Return list of Fastq files in the set
        """
        return list(filter(lambda fq: fq is not None,
                           self._fastqs))

class QCReport(Document):
    """
    Create a QC report document for one or more projects

    Example usage:

    >>> report = QCReport(project)
    >>> report.write("qc_report.html")

    To control the fields written to the summary table, specify
    a list of field names via the 'summary_fields' argument.
    Valid field names are:

    - sample: sample name
    - fastq: Fastq name
    - fastqs: Fastq R1/R2 names
    - reads: number of reads
    - read_lengths: length of reads
    - fastqc_[read]: FastQC mini-plot for [read] (r1,r2,...)
    - boxplot_[read]: FastQC per-base-quality mini-boxplot' for [read]
    - screens_[read]: FastQScreen mini-plots for [read]
    - strandedness: 'forward', 'reverse' or 'unstranded' for pair
    - cellranger_count: 'cellranger count' outputs for each sample

    To control the elements written to the reports for each Fastq
    pair, specify a list of element names via the 'report_attrs'
    argument. Valid element names are:

    - fastqc: FastQC report
    - fastq_screen: FastQCScreen report
    - program_versions: program versions
    """
    # Field descriptions for summary table
    field_descriptions = {
        'sample': 'Sample',
        'fastq' : 'Fastq',
        'fastqs': 'Fastqs',
        'reads': '#reads',
        'read_lengths': 'Length',
        'fastqc_r1': 'FastQC[R1]',
        'boxplot_r1': 'Boxplot[R1]',
        'screens_r1': 'Screens[R1]',
        'fastqc_r2': 'FastQC[R2]',
        'boxplot_r2': 'Boxplot[R2]',
        'screens_r2': 'Screens[R2]',
        'fastqc_r3': 'FastQC[R3]',
        'boxplot_r3': 'Boxplot[R3]',
        'screens_r3': 'Screens[R3]',
        'strandedness': 'Strand',
        'cellranger_count': 'Cellranger count',
    }
    # Titles for metadata items
    metadata_titles = {
        'project_id': 'Project ID',
        'run_id': 'Run ID',
        'run': 'Run name',
        'user': 'User',
        'PI': 'PI',
        'library_type': 'Library type',
        'sequencer_model': 'Sequencer model',
        'single_cell_platform': 'Single cell preparation platform',
        'number_of_cells': 'Number of cells',
        'organism': 'Organism',
        'protocol': 'QC protocol',
        'cellranger_reference': 'Cellranger reference dataset',
        'multiqc': 'MultiQC report',
        'icell8_stats': 'ICELL8 statistics',
        'icell8_report': 'ICELL8 processing report',
    }
    # Software packages and names
    software_packages = ['bcl2fastq',
                         'cellranger',
                         'cellranger-atac',
                         'cellranger-arc',
                         'spaceranger',
                         'fastqc',
                         'fastq_screen',
                         'fastq_strand',]
    software_names = {
        'bcl2fastq': 'Bcl2fastq',
        'cellranger': 'Cellranger',
        'cellranger-atac': 'Cellranger ATAC',
        'cellranger-arc': 'Cellranger ARC',
        'spaceranger': 'Spaceranger',
        'fastqc': 'FastQC',
        'fastq_screen': 'FastqScreen',
        'fastq_strand': 'FastqStrand',
    }
    def __init__(self,projects,title=None,qc_dir=None,report_attrs=None,
                 summary_fields=None,relpath=None,data_dir=None):
        """
        Create a new QCReport instance

        Arguments:
          projects (AnalysisProject): list of projects to
             report QC for
          title (str): title for the report (defaults to
            "QC report: <PROJECT_NAME>")
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
          report_attrs (list): list of elements to report for
            each Fastq pair
          summary_fields (list): list of fields to report for
            each sample in the summary table
          relpath (str): if set then make link paths
            relative to 'relpath'
          data_dir (str): if set then copy external data
            files to this directory and make link paths
            to these copies; relative path will be treated
            as a subdirectory of the project
        """
        # Convert projects to QCProjects
        projects = [QCProject(p,qc_dir=qc_dir) for p in projects]
        # Flag to indicate if report references multiple projects
        self.multi_project = (len(projects) > 1)
        # Flag to indicate if any project has single cell data
        self.has_single_cell = any([p.is_single_cell for p in projects])
        # Primary project
        primary_project = projects[0]
        # Set up title
        if title is None:
            if self.multi_project:
                title = "QC report: %s &amp; %s" % (
                    ', '.join([p.id for p in projects[0:-1]]),
                    projects[-1].id)
            else:
                title = "QC report: %s" % primary_project.name
        # Initialise superclass
        Document.__init__(self,title)
        # Collect initial QC data across all projects
        self.output_files = []
        self.outputs = set()
        for project in projects:
            # Output files
            if project.output_files:
                self.output_files.extend([f for f in project.output_files])
            # Outputs
            for output in project.outputs:
                self.outputs.add(output)
        self.outputs = sorted(list(self.outputs))
        # Software packages
        self.software = []
        for pkg in self.software_packages:
            for project in projects:
                if project.software_info(pkg):
                    self.software.append(pkg)
                    break
        # Use data directory?
        if data_dir is not None:
            if not os.path.abspath(data_dir):
                data_dir = os.path.join(primary_project.dirn,data_dir)
            data_dir = os.path.normpath(data_dir)
        self.data_dir = data_dir
        # Relative paths
        if relpath is not None:
            relpath = os.path.normpath(os.path.abspath(relpath))
        self.relpath = relpath
        # Status of report
        self.status = True
        # Initialise tables
        self._init_metadata_table(projects)
        self._init_software_table()
        # Initialise report sections
        self.preamble = self._init_preamble_section()
        self.warnings = self._init_warnings_section()
        self.summary = self._init_summary_section(project)
        # Initialise data directory
        self._init_data_dir(projects)
        # Report each project
        for project in projects:
            # Report QC data for project
            print("Project: %s" % project.name)
            if project.samples:
                print("Samples found:")
                for sample in project.samples:
                    print("\t- %s" % sample)
            else:
                logger.warning("%s: no samples found" % project.name)
            if project.fastqs:
                print("Fastqs referenced:")
                for fastq in project.fastqs:
                    print("\t- %s" % fastq)
                print("Reads found:")
                for read in project.reads:
                    print("\t- %s" % read)
            else:
                logger.warning("%s: no Fastqs referenced" % project.name)
            if project.outputs:
                print("Available QC outputs:")
                for output in project.outputs:
                    print("\t- %s" % output)
            else:
                logger.warning("%s: no QC outputs found" % project.name)
            if project.software:
                print("Software versions:")
                for package in project.software:
                    print("\t- %s: %s" %
                          (package,
                           ','.join(project.software[package])))
            # Fields to report in summary table
            if not summary_fields:
                if len(project.reads) > 1:
                    summary_fields_ = ['sample',
                                       'fastqs',
                                       'reads',
                                       'read_lengths',]
                else:
                    summary_fields_ = ['sample',
                                       'fastq',
                                       'reads',
                                       'read_lengths',]
                if 'strandedness' in project.outputs:
                    summary_fields_.append('strandedness')
                for read in project.reads:
                    if ('fastqc_%s' % read) in project.outputs:
                        summary_fields_.append('fastqc_%s' % read)
                        summary_fields_.append('boxplot_%s' % read)
                    if ('screens_%s' % read) in project.outputs:
                        summary_fields_.append('screens_%s' % read)
                if 'cellranger_count' in project.outputs:
                    summary_fields_.append('cellranger_count')
            # Attributes to report for each sample
            if report_attrs is None:
                report_attrs_ = ['fastqc',
                                'fastq_screen',]
                if 'strandedness' in project.outputs:
                    report_attrs_.append('strandedness')
            # Add data for this project to the report
            print("Adding project '%s' to the report..." % project.name)
            self.report_metadata(project)
            self.report_software(project)
            self.report_comments(project)
            # Create a new summary table
            summary_table = self.add_summary_table(project,summary_fields_)
            # Report each sample
            for sample in project.samples:
                self.report_sample(project,sample,report_attrs_,
                                   summary_table,summary_fields_)
        # Report the status
        self.report_status()

    def _init_metadata_table(self,projects):
        """
        Internal: set up a table for project metadata

        Associated CSS class is 'metadata'
        """
        # Identify metadata items
        metadata_items = ['run_id',
                          'run',
                          'sequencer_model',
                          'user',
                          'PI',
                          'library_type',
                          'organism',
                          'protocol',]
        if self.has_single_cell:
            for item in ('single_cell_platform',
                         'number_of_cells',):
                metadata_items.insert(metadata_items.index('organism'),
                                      item)
        if 'cellranger_count' in self.outputs:
            metadata_items.append('cellranger_reference')
        if 'multiqc' in self.outputs:
            metadata_items.append('multiqc')
        if 'icell8_stats' in self.outputs:
            metadata_items.append('icell8_stats')
        if 'icell8_report' in self.outputs:
            metadata_items.append('icell8_report')
        if self.multi_project:
            metadata_items[metadata_items.index('run_id')] = 'project_id'
        # Make table with one column per project
        columns = ['item']
        for project in projects:
            columns.append(project.id)
        metadata_table = Table(columns)
        metadata_table.no_header()
        metadata_table.add_css_classes('metadata')
        # Add rows for metadata items
        for item in metadata_items:
            metadata_table.add_row(item=self.metadata_titles[item])
        # Store table and metadata items as attribute
        self.metadata_table = metadata_table
        self.metadata_items = metadata_items

    def _init_software_table(self):
        """
        Internal: set up a table for software information

        Associated CSS class is 'metadata'
        """
        # Determine what software packages were used across
        # all projects
        software_table = Table(('program','version',))
        software_table.no_header()
        software_table.add_css_classes('metadata')
        for pkg in self.software:
            software_table.add_row(
                program=self.software_names[pkg],
                version=None)
        self.software_table = software_table

    def _init_preamble_section(self):
        """
        Internal: set up a "preamble" section

        Associated name is 'preamble'
        """
        preamble = self.add_section(name='preamble')
        preamble.add("Report generated by auto_process %s on %s" %
                     (get_version(),time.asctime()))
        return preamble

    def _init_summary_section(self,project):
        """
        Internal: set up a summary section for the report

        Associated name is 'summary'
        """
        summary = self.add_section("Summary",name='summary',
                                   css_classes=("summary",))
        # Add subsections inside a container
        info = summary.add_subsection()
        general_info = info.add_subsection("General information",
                                           css_classes=("info",))
        general_info.add(self.metadata_table)
        software_info = info.add_subsection("Software",
                                           css_classes=("info",))
        software_info.add(self.software_table)
        # Add an empty section to clear HTML floats
        clear = summary.add_subsection(css_classes=("clear",))
        # Add additional subsections for comments etc
        info = summary.add_subsection()
        # Add comments section
        self.comments = info.add_subsection("Comments",
                                            css_classes=("info",))
        # Add an empty section to clear HTML floats
        clear = summary.add_subsection(css_classes=("clear",))
        return summary

    def _init_warnings_section(self):
        """
        Internal: creates a section which is displayed if there
        are warnings
        """
        warnings = self.add_section(name='warnings',
                                    css_classes=("warnings",))
        warnings.add(Para(WarningIcon(size=50),
                          "There are issues with this project"))
        return warnings

    def _init_data_dir(self,projects):
        """
        Internal: initialise the data directory

        Creates directory, copies QC artefacts and updates the
        list of output files associated with the report
        """
        if not self.data_dir:
            return
        # Remove pre-existing data directory
        if os.path.exists(self.data_dir):
            print("Removing existing data directory: %s" %
                  self.data_dir)
            shutil.rmtree(self.data_dir)
        # (Re)create directory
        print("Creating data directory: %s" % self.data_dir)
        os.mkdir(self.data_dir)
        # Copy QC artefacts from each project
        output_files = []
        for project in projects:
            # Set up a project-specific top-level directory
            project_data_dir = os.path.join(self.data_dir,
                                            sanitize_name(project.id))
            if not os.path.exists(project_data_dir):
                logging.debug("Creating subdirectory %s" %
                              project_data_dir)
                os.makedirs(project_data_dir)
            for output_file in project.output_files:
                # Convert file paths to be relative to project dir
                output_file = os.path.relpath(output_file,
                                              project.dirn)
                # Make a copy of the file in the data dir
                subdir = os.path.dirname(output_file)
                # File includes a leading directory path so
                # build a copy of this first
                if subdir:
                    subdir = os.path.join(project_data_dir,subdir)
                    if not os.path.exists(subdir):
                        os.makedirs(subdir)
                # Copy the file
                src = os.path.join(project.dirn,output_file)
                dest = os.path.join(project_data_dir,output_file)
                logging.debug("-- copying %s to %s" % (output_file,
                                                       dest))
                shutil.copyfile(src,dest)
                # Add to the updated list of output files
                output_files.append(dest)
        # Update the list of files
        self.output_files = output_files

    def add_summary_table(self,project,fields):
        """
        Create a new table for summarising samples from a project

        Associated CSS classes are 'summary' and 'fastq_summary'
        """
        # Add header For multi-project
        if self.multi_project:
            summary = self.summary.add_subsection(
                project.id,
                name=sanitize_name(project.id))
        else:
            summary = self.summary
        # Create the table
        summary_tbl = Table(fields,**self.field_descriptions)
        summary_tbl.add_css_classes('summary','fastq_summary')
        # Append to the summary section
        summary.add(
            "%d sample%s | %d fastq%s" % (
                len(project.samples),
                ('s' if len(project.samples) != 1 else ''),
                len(project.fastqs),
                ('s' if len(project.fastqs) != 1 else ''))
        )
        summary.add(summary_tbl)
        return summary_tbl

    def report_metadata(self,project):
        """
        Report the project metadata

        Adds entries for the project metadata to the "metadata"
        table in the report
        """
        # Determine the root directory for QC outputs
        if self.data_dir:
            project_data_dir = os.path.join(self.data_dir,
                                            sanitize_name(project.id))
        else:
            project_data_dir = project.dirn
        if self.relpath:
            project_data_dir = os.path.relpath(project_data_dir,
                                               self.relpath)
        # Add metadata items
        for idx,item in enumerate(self.metadata_items):
            # Try to acquire the value from QC metadata
            try:
                value = project.qc_info[item]
            except KeyError:
                # Fall back to project metadata
                try:
                    if project.info[item]:
                        value = project.info[item]
                    else:
                        # No value set, skip this item
                        continue
                except KeyError:
                    # Additional non-metadata items, or items
                    # requiring additional processing
                    if item == 'project_id':
                        value = Link(project.id,
                                     "#%s" % sanitize_name(project.id))
                    elif item == 'run_id':
                        value = project.run_id
                        if value is None:
                            # Unable to determine run id
                            continue
                    elif item == 'cellranger_reference':
                        path = project.qc_info.cellranger_refdata
                        if path is None:
                            # No reference dataset
                            continue
                        if os.path.dirname(path):
                            value = "...%s%s" % (os.sep,
                                                 os.path.basename(path))
                        else:
                            value = path
                    elif item == 'multiqc':
                        multiqc_report = "multi%s_report.html" \
                                         % os.path.basename(project.qc_dir)
                        value = Link(multiqc_report,
                                     os.path.join(project_data_dir,
                                                  multiqc_report))
                    elif item == 'icell8_stats':
                        value = Link("icell8_stats.xlsx",
                                     os.path.join(project_data_dir,
                                                  "stats",
                                                  "icell8_stats.xlsx"))
                    elif item == 'icell8_report':
                        value = Link("icell8_processing.html",
                                     os.path.join(project_data_dir,
                                                  "icell8_processing.html"))
                    else:
                        raise Exception("Unrecognised item to report: '%s'"
                                        % item)
            # Update the value in the metadata table
            self.metadata_table.set_value(idx,
                                          project.id,
                                          value)

    def report_software(self,project):
        """
        Report the software versions used in the processing & QC

        Adds entries for the software versions to the "software"
        table in the report
        """
        # Report software packages
        for pkg in self.software_packages:
            # Acquire the value
            value = project.software_info(pkg)
            if value is None:
                continue
            # Get row index in the metadata table
            idx = self.software.index(pkg)
            self.software_table.set_value(idx,
                                          "version",
                                          value)

    def report_comments(self,project):
        """
        Report the comments associated with a project

        Adds the comments from the project metadata as
        a list to the comments section.

        Arguments:
          project (QCProject): project to report
        """
        comments_list = List()
        try:
            if project.comments:
                for comment in project.comments.split(';'):
                    comments_list.add_item(comment.strip())
            else:
                # Drop out with exception
                raise AttributeError
        except AttributeError:
            comments_list.add_item("N/A")
        if self.multi_project:
            self.comments.add("<p><b>%s</b></p>" % project.id)
        self.comments.add(comments_list)

    def report_sample(self,project,sample,report_attrs,summary_table,
                      summary_fields):
        """
        Report the QC for a sample

        Reports the QC for the sample and Fastqs to
        the summary table and appends a section with
        detailed reports to the document.

        Arguments:
          project (QCProject): project to report
          sample (str): name of sample to report
          report_attrs (list): list of elements to report for
            each set of Fastqs
          summary_table (Table): summary table to report
            each sample in
          summary_fields (list): list of fields to report
            for each sample in the summary table
        """
        # Create a unique name and title
        if self.multi_project:
            sample_name = "sample_%s_%s" % (sanitize_name(project.id),
                                            sample)
            sample_title = "Sample: %s/%s" % (project.id,sample)
        else:
            sample_name = "sample_%s" % sample
            sample_title = "Sample: %s" % sample
        # Create a new section
        sample_report = self.add_section(
            sample_title,
            name=sample_name,
            css_classes=('sample',))
        # Get Fastq groups
        fastqs = sorted(list(
            filter(lambda fq:
                   project.fastq_attrs(fq).sample_name == sample,
                   project.fastqs)))
        fastq_groups = group_fastqs_by_name(fastqs,project.fastq_attrs)
        # Location of QC artefacts
        qc_dir = project.qc_dir
        if self.data_dir:
            qc_dir = os.path.join(self.data_dir,
                                  sanitize_name(project.id),
                                  os.path.basename(qc_dir))
        # Report number of fastqs and reads
        reads = QCReportFastqGroup(fastq_groups[0],
                                   qc_dir=qc_dir,
                                   project_id=project.id,
                                   fastq_attrs=project.fastq_attrs).reads
        if len(reads) == 1:
            sample_report.add("%d %s Fastq%s" %
                              (len(fastq_groups),
                               reads[0].upper(),
                               's' if len(fastq_groups) > 1 else ''))
        elif len(reads) == 2:
            sample_report.add("%d %s Fastq pair%s" %
                              (len(fastq_groups),
                               '/'.join([r.upper() for r in reads]),
                              's' if len(fastq_groups) > 1 else ''))
        else:
            sample_report.add("%d %s Fastq group%s" %
                              (len(fastq_groups),
                               '/'.join([r.upper() for r in reads]),
                              's' if len(fastq_groups) > 1 else ''))
        # Keep track of the first line in the summary
        # table, as per-sample metrics (and name)
        # should only be reported on the first line
        first_line = True
        # Report each Fastq group
        for fastqs in fastq_groups:
            # Report Fastq pair
            fastq_group = QCReportFastqGroup(fastqs,
                                             qc_dir=qc_dir,
                                             project_id=project.id,
                                             fastq_attrs=project.fastq_attrs)
            fastq_group.report(sample_report,
                               attrs=report_attrs,
                               relpath=self.relpath)
            # Add line in summary table
            if sample is not None:
                if first_line:
                    # Only display sample name on first line
                    idx = summary_table.add_row(
                        sample=Link(sample,
                                    sample_report))
                else:
                    # Don't display sample name for subsequent
                    # lines in summary table
                    idx = summary_table.add_row(sample="&nbsp;")
            else:
                idx = summary_table.add_row(sample="&nbsp;")
            status = fastq_group.update_summary_table(
                summary_table,idx=idx,
                fields=summary_fields,
                relpath=self.relpath,
                skip_sample_metrics=(not first_line))
            if not status:
                # Update flag to indicate problems with the
                # report
                self.status = False
            # Update flag to indicate no longer on
            # first line for this sample
            first_line = False

    def report_status(self):
        """
        Set the visibility of the "warnings" section
        """
        if self.status:
            # Turn off display of warnings section
            self.warnings.add_css_classes("hide")

class QCReportFastqGroup(object):
    """
    Utility class for reporting the QC for a Fastq group

    Provides the following properties:

    reads: list of read ids e.g. ['r1','r2']
    fastqs: dictionary mapping read ids to Fastq paths
    reporters: dictionary mapping read ids to QCReportFastq
      instances
    paired_end: whether FastqGroup is paired end
    fastq_strand_txt: location of associated Fastq_strand
      output

    Provides the following methods:

    strandedness: fetch strandedness data for this group
    ustrandplot: return mini-strand stats summary plot
    report_strandedness: write report for strandedness
    report: write report for the group
    update_summary_table: add line to summary table for
      the group
    """
    def __init__(self,fastqs,qc_dir,project_id=None,
                 fastq_attrs=AnalysisFastq):
        """
        Create a new QCReportFastqGroup

        Arguments:
          fastqs (list): list of paths to Fastqs in the
            group
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
          project_id (str): identifier for the project
          fastq_attrs (BaseFastqAttrs): class for extracting
            data from Fastq names
        """
        self.qc_dir = qc_dir
        self.project_id = project_id
        self.fastq_attrs = fastq_attrs
        # Assign fastqs to reads
        self.fastqs = defaultdict(lambda: None)
        self.reporters = defaultdict(lambda: None)
        self.reads = set()
        for fastq in fastqs:
            fq = self.fastq_attrs(fastq)
            if fq.is_index_read:
                read = 'i'
            else:
                read = 'r'
            read = "%s%d" % (read,
                             fq.read_number
                             if fq.read_number is not None else 1)
            self.fastqs[read] = fastq
            self.reporters[read] = QCReportFastq(fastq,
                                                 self.qc_dir,
                                                 self.project_id,
                                                 fastq_attrs=
                                                 self.fastq_attrs)
            self.reads.add(read)
        self.reads = sorted(list(self.reads))

    @property
    def paired_end(self):
        """
        True if pair consists of R1/R2 files
        """
        return ('r1' in self.fastqs and 'r2' in self.fastqs)

    @property
    def fastq_strand_txt(self):
        """
        Locate output from fastq_strand (None if not found)
        """
        fastq_strand_txt = None
        for fq in (self.fastqs['r1'],self.fastqs['r2']):
            if fq is not None:
                fastq_strand_txt = os.path.join(
                    self.qc_dir,
                    fastq_strand_output(fq))
                if os.path.isfile(fastq_strand_txt):
                    return fastq_strand_txt
        return None

    def strandedness(self):
        """
        Return strandedness from fastq_strand.py
        """
        # Locate strandedness data
        txt = self.fastq_strand_txt
        # No file found
        if txt is None:
            return "!!!No strandedness data!!!"
        # Return stats
        strandedness = Fastqstrand(txt)
        output = []
        for genome in strandedness.genomes:
            try:
                ratio = "%.2f" % strandedness.stats[genome].ratio
            except TypeError:
                ratio = "NaN"
            output.append("%s: F%.2f%% | R%.2f%% (%s)" %
                          (genome,
                           strandedness.stats[genome].forward,
                           strandedness.stats[genome].reverse,
                           ratio))
        return "\n".join(output)

    def ustrandplot(self):
        """
        Return a mini-strand stats summary plot
        """
        txt = self.fastq_strand_txt
        return ustrandplot(txt,inline=True,dynamic=True)

    def report_strandedness(self,document):
        """
        Report the strandedness outputs to a document

        Creates a new subsection called "Strandedness"
        with a table of the strandedness determination
        outputs.

        Arguments:
          document (Section): section to add report to
        """
        strandedness_report = document.add_subsection(
            "Strandedness",
            name="strandedness_%s" % self.reporters['r1'].safe_name)
        strandedness_report.add_css_classes("strandedness")
        # Locate strandedness
        txt = self.fastq_strand_txt
        # No file found
        if txt is None:
            strandedness_report.add(WarningIcon(),
                                    "No strandedness data available")
        else:
            strandedness = Fastqstrand(txt)
            strandedness_report.add("Strandedness statistics from "
                                    "<tt>fastq_strand %s</tt>:"
                                    % strandedness.version)
            strandedness_tbl = Table(("genome",
                                      "forward",
                                      "reverse",
                                      "ratio",
                                      "strand"),
                                     genome="Genome",
                                     forward="Forward %",
                                     reverse="Reverse %",
                                     ratio="Ratio (forward/reverse)",
                                     strand="Strandedness *")
            strandedness_tbl.add_css_classes("summary")
            for genome in strandedness.genomes:
                try:
                    ratio = "%.2f" % strandedness.stats[genome].ratio
                except TypeError:
                    ratio = "NaN"
                strandedness_tbl.add_row(
                    genome=genome,
                    forward=strandedness.stats[genome].forward,
                    reverse=strandedness.stats[genome].reverse,
                    ratio=ratio,
                    strand=strandedness.stats[genome].strandedness)
            strandedness_report.add(strandedness_tbl)
            strandedness_report.add("* Strandedness is 'reverse' if "
                                    "ratio &lt;0.2, 'forward' if ratio "
                                    "is &gt;5, 'unstranded' otherwise "
                                    "(or undetermined if 'NaN')")
        return strandedness_report

    def report(self,sample_report,attrs=None,relpath=None):
        """
        Add report for Fastq group to a document section

        Creates a new subsection in 'sample_report' for the
        Fastq pair, within which are additional subsections
        for each Fastq file.

        The following 'attributes' can be reported for each
        Fastq:

        - fastqc
        - fastq_screen
        - program versions

        By default all attributes are reported.

        Arguments:
          sample_report (Section): section to add the report
            to
          attrs (list): optional list of custom 'attributes'
            to report
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        # Attributes to report
        if attrs is None:
            attrs = ('fastqc',
                     'fastq_screen',
                     'program_versions',
                     'strandedness')
        # Add container section for Fastq pair
        fastqs_report = sample_report.add_subsection(css_classes=('fastqs',))
        # Create sections for individual Fastqs
        for read in self.reads:
            fq = self.reporters[read]
            if fq is None:
                continue
            fq_report = fastqs_report.add_subsection(fq.name,
                                                     name=fq.safe_name,
                                                     css_classes=('fastq',))
            # Add reports for each requested 'attribute'
            for attr in attrs:
                if attr == "fastqc":
                    # FastQC outputs
                    fq.report_fastqc(fq_report,relpath=relpath)
                elif attr == "fastq_screen":
                    # FastQScreens
                    fq.report_fastq_screens(fq_report,relpath=relpath)
                elif attr == "program_versions":
                    # Versions of programs used
                    fq.report_program_versions(fq_report)
                elif attr == "strandedness":
                    # Strandedness - handle separately
                    pass
                else:
                    raise KeyError("'%s': unrecognised reporting element "
                                   % attr)
        # Sections for pairwise data
        if "strandedness" in attrs:
            # Strandedness
            self.report_strandedness(fastqs_report)
        # Add an empty section to clear HTML floats
        clear = fastqs_report.add_subsection(css_classes=("clear",))

    def update_summary_table(self,summary_table,idx=None,fields=None,
                             relpath=None,skip_sample_metrics=False):
        """
        Add a line to a summary table reporting a Fastq group

        Creates a new line in 'summary_table' (or updates an
        existing line) for the Fastq pair, adding content for
        each specified field.

        The following fields can be reported for each Fastq
        pair:

        - fastqs (if paired-end)
        - fastq (if single-end)
        - reads
        - read_lengths
        - boxplot_r1
        - boxplot_r2
        - boxplot_r3
        - fastqc_r1
        - fastqc_r2
        - fastqc_r3
        - screens_r1
        - screens_r2
        - screens_r3
        - strandedness
        - cellranger_count

        Arguments:
          summary_table (Table): table to add the summary to
          idx (int): if supplied then indicates which existing
            table row to update (if None then a new row is
            appended)
          fields (list): list of custom fields to report
          relpath (str): if set then make link paths
            relative to 'relpath'
          skip_sample_metrics (bool): if True then don't report
            values for 'sample-level' metrics (e.g. cellranger
            count outputs)

        Returns:
          Boolean: True if report didn't contain any issues,
            False otherwise.
        """
        # Flag indicating issues
        has_problems = False
        # Fields to report
        if fields is None:
            if self.paired_end:
                fields = ('fastqs',
                          'reads',
                          'read_lengths',
                          'boxplot_r1','boxplot_r2',
                          'fastqc_r1','fastqc_r2',
                          'screens_r1','screens_r2')
            else:
                fields = ('fastq',
                          'reads',
                          'read_lengths',
                          'boxplot_r1',
                          'fastqc_r1',
                          'screens_r1')
        else:
            # Drop fields for reads that aren't present
            # For example if there are a mixture of
            # single and paired-end Fastqs
            updated_fields = []
            for field in fields:
                if field[:-1].endswith("_r"):
                    read = field.split("_")[-1]
                    if read not in self.reads:
                        logging.warning("Dropping '%s' from summary table"
                                        % field)
                        continue
                updated_fields.append(field)
            fields = updated_fields
        # Add row to summary table
        if idx is None:
            idx = summary_table.add_row()
        # Populate with data
        for field in fields:
            try:
                if field == "sample":
                    logger.debug("'sample' ignored")
                    continue
                elif field == "fastq" or field == "fastqs":
                    value = []
                    for read in self.reads:
                        value.append(Link(self.reporters[read].name,
                                          "#%s" % self.reporters[read].safe_name))
                    value = "<br />".join([str(x) for x in value])
                elif field == "reads":
                    value = pretty_print_reads(
                        self.reporters[self.reads[0]].fastqc.data.basic_statistics(
                            'Total Sequences'))
                elif field == "read_lengths":
                    value = []
                    for read in self.reads:
                        value.append(Link(
                            self.reporters[read].fastqc.data.basic_statistics(
                                'Sequence length'),
                            self.reporters[read].fastqc.summary.link_to_module(
                                'Sequence Length Distribution',
                                relpath=relpath)))
                    value = "<br />".join([str(x) for x in value])
                elif field.startswith("boxplot_"):
                    read = field.split('_')[-1]
                    value = Img(self.reporters[read].uboxplot(),
                                href="#boxplot_%s" %
                                self.reporters[read].safe_name)
                elif field.startswith("fastqc_"):
                    read = field.split('_')[-1]
                    value = Img(self.reporters[read].ufastqcplot(),
                                href="#fastqc_%s" %
                                self.reporters[read].safe_name,
                                title=self.reporters[read].fastqc_summary())
                elif field.startswith("screens_"):
                    read = field.split('_')[-1]
                    value = Img(self.reporters[read].uscreenplot(),
                                href="#fastq_screens_%s" %
                                self.reporters[read].safe_name)
                elif field == "strandedness":
                    value = Img(self.ustrandplot(),
                                href="#strandedness_%s" %
                                self.reporters[self.reads[0]].safe_name,
                                title=self.strandedness())
                elif field == "cellranger_count":
                    if skip_sample_metrics:
                        value = "&nbsp;"
                    else:
                        web_summary = self.reporters[read].\
                                      cellranger_count.web_summary
                        if relpath:
                            web_summary = os.path.relpath(web_summary,
                                                          relpath)
                        value = Link(
                            self.reporters[read].cellranger_count.sample_name,
                            web_summary)
                else:
                    raise KeyError("'%s': unrecognised field for summary "
                                   "table" % field)
                # Put value into the table
                summary_table.set_value(idx,field,value)
            except KeyError as ex:
                # Assume this is an unrecognised field
                logger.error("%s" % ex)
                raise ex
            except Exception as ex:
                # Encountered an exception trying to acquire the value
                # for the field
                fastqs = ", ".join([self.reporters[r].name
                                    for r in self.reads])
                logger.warning("Exception setting '%s' in summary table "
                              "for Fastq group { %s }: %s"
                               % (field,
                                  fastqs,
                                  ex))
                # Put error value into the table
                summary_table.set_value(idx,field,
                                        WarningIcon("Unable to get value for "
                                                    "%s for %s" % (field,
                                                                   fastqs),
                                                     size=25))
                # Update flag
                has_problems = True
        # Return the status
        return (not has_problems)

class QCReportFastq(object):
    """
    Provides interface to QC outputs for Fastq file

    Provides the following attributes:

    name: basename of the Fastq
    path: path to the Fastq
    safe_name: name suitable for use in HTML links etc
    sample_name: sample name derived from the Fastq basename
    fastqc: Fastqc instance
    fastq_screen.names: list of FastQScreen names
    fastq_screen.SCREEN.description: description of SCREEN
    fastq_screen.SCREEN.png: associated PNG file for SCREEN
    fastq_screen.SCREEN.txt: associated TXT file for SCREEN
    fastq_screen.SCREEN.version: associated version for SCREEN
    program_versions.NAME: version of package NAME
    cellranger_count.web_summary: path to cellranger count web_summary.html
    cellranger_count.metrics_cvs: path to cellranger count metrics.csv

    Provides the following methods:

    report_fastqc
    report_fastq_screens
    report_program_versions
    uboxplot
    ufastqcplot
    uscreenplot
    """
    def __init__(self,fastq,qc_dir,project_id=None,fastq_attrs=AnalysisFastq):
        """
        Create a new QCReportFastq instance

        Arguments:
          fastq (str): path to Fastq file
          qc_dir (str): path to QC directory
          project_id (str): identifier for the parent project
          fastq_attrs (BaseFastqAttrs): class for extracting
            data from Fastq names
        """
        # Source data
        logging.debug("******** QCREPORTFASTQ ************")
        self.name = os.path.basename(fastq)
        self.path = os.path.abspath(fastq)
        self.project_id = project_id
        self.safe_name = strip_ngs_extensions(self.name)
        if self.project_id:
            self.safe_name = "%s_%s" % (sanitize_name(self.project_id),
                                        self.safe_name)
        self.fastq_attrs = fastq_attrs
        # Sample name
        self.sample_name = fastq_attrs(self.name).sample_name
        logging.debug("Name  : %s" % self.name)
        logging.debug("Path  : %s" % self.path)
        logging.debug("ID    : %s" % self.project_id)
        logging.debug("Nsafe : %s" % self.safe_name)
        logging.debug("QCDir : %s" % qc_dir)
        # FastQC
        logging.debug("Fastqc: %s" %
                      os.path.join(qc_dir,fastqc_output(fastq)[0]))
        try:
            self.fastqc = Fastqc(os.path.join(
                qc_dir,fastqc_output(fastq)[0]))
        except Exception as ex:
            self.fastqc = None
        # Fastqscreen
        self.fastq_screen = AttributeDictionary()
        self.fastq_screen['names'] = list()
        for name in FASTQ_SCREENS:
            png,txt = fastq_screen_output(fastq,name)
            png = os.path.join(qc_dir,png)
            txt = os.path.join(qc_dir,txt)
            if os.path.exists(png) and os.path.exists(txt):
                self.fastq_screen['names'].append(name)
                self.fastq_screen[name] = AttributeDictionary()
                self.fastq_screen[name]["description"] = \
                                        name.replace('_',' ').title()
                self.fastq_screen[name]["png"] = png
                self.fastq_screen[name]["txt"] = txt
                self.fastq_screen[name]["version"] = Fastqscreen(txt).version
        # Program versions
        self.program_versions = AttributeDictionary()
        if self.fastqc is not None:
            fastqc_version = self.fastqc.version
        else:
            fastqc_version = WarningIcon(size=20)
        self.program_versions['fastqc'] = fastqc_version
        fastq_screen_versions = list(
            set([self.fastq_screen[s].version
                 for s in self.fastq_screen.names]))
        if fastq_screen_versions:
            fastq_screen_versions = ','.join(sorted(fastq_screen_versions))
        else:
            fastq_screen_versions = WarningIcon(size=20)
        self.program_versions['fastq_screen'] = fastq_screen_versions
        # Cellranger count outputs
        self.cellranger_count = CellrangerCount(
            os.path.join(qc_dir,
                         "cellranger_count",
                         self.sample_name))

    def report_fastqc(self,document,relpath=None):
        """
        Report the FastQC outputs to a document

        Creates a new subsection called "FastQC" with
        a copy of the FastQC sequence quality boxplot and
        a summary table of the results from each FastQC
        module.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        fastqc_report = document.add_subsection("FastQC")
        if self.fastqc:
            # FastQC quality boxplot
            fastqc_report.add("Per base sequence quality boxplot:")
            boxplot = Img(self.fastqc.quality_boxplot(inline=True),
                          height=250,
                          width=480,
                          href=self.fastqc.summary.link_to_module(
                              'Per base sequence quality',
                              relpath=relpath),
                          name="boxplot_%s" % self.safe_name)
            fastqc_report.add(boxplot)
            # FastQC summary table
            fastqc_report.add("FastQC summary:")
            fastqc_tbl = Target("fastqc_%s" % self.safe_name)
            fastqc_report.add(fastqc_tbl,
                              self.fastqc.summary.html_table(relpath=relpath))
            if relpath:
                fastqc_html_report = os.path.relpath(self.fastqc.html_report,
                                                     relpath)
            else:
                fastqc_html_report = self.fastqc.html_report
            fastqc_report.add("%s for %s" % (Link("Full FastQC report",
                                                  fastqc_html_report),
                                             self.name))
        else:
            fastqc_report.add(WarningIcon(),"No FastQC data available")
        return fastqc_report

    def report_fastq_screens(self,document,relpath=None):
        """
        Report the FastQScreen outputs to a document

        Creates a new subsection called "Screens" with
        copies of the screen plots for each screen and
        links to the "raw" text files.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        screens_report = document.add_subsection("Screens",
                                                 name="fastq_screens_%s" %
                                                 self.safe_name)
        if not self.fastq_screen.names:
            screens_report.add(WarningIcon(),"No screens found")
            return screens_report
        raw_data = list()
        for name in self.fastq_screen.names:
            # Title for the screen
            screens_report.add(self.fastq_screen[name].description)
            # Gather data for screen
            png = self.fastq_screen[name].png
            txt = self.fastq_screen[name].txt
            if relpath:
                png_href = os.path.relpath(png,relpath)
                txt_href = os.path.relpath(txt,relpath)
            else:
                png_href = png
                txt_href = txt
            # Screen plot (PNG)
            if os.path.exists(png):
                screens_report.add(Img(encode_png(png),
                                       height=250,
                                       href=png_href))
            else:
                screens_report.add("!!!No FastqScreen plot available!!!")
            # Link to raw data (TXT)
            if os.path.exists(txt):
                raw_data.append(
                    Link(self.fastq_screen[name].description,
                         txt_href).html())
            else:
                screens_report.add("!!!No FastqScreen data available!!!")
        # Add links to raw data (TXT files)
        screens_report.add("Raw screen data: " +
                           "| ".join(raw_data))
        return screens_report

    def report_program_versions(self,document):
        """
        Report the program versions to a document

        Creates a new subsection called "Program versions"
        with a table listing the versions of the QC
        programs.

        Arguments:
          document (Section): section to add report to
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        versions = document.add_subsection("Program versions")
        programs = Table(("Program","Version"))
        programs.add_css_classes("programs","summary")
        programs.add_row(Program='fastqc',
                         Version=self.program_versions.fastqc)
        programs.add_row(Program='fastq_screen',
                         Version=self.program_versions.fastq_screen)
        versions.add(programs)
        return versions

    def fastqc_summary(self):
        """
        Return plaintext version of the FastQC summary
        """
        output = []
        for m in self.fastqc.summary.modules:
            output.append("%s: %s" %
                          (self.fastqc.summary.status(m),m))
        return "\n".join(output)

    def uboxplot(self,inline=True):
        """
        Return a mini-sequence quality boxplot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return uboxplot(self.fastqc.data.path,inline=inline)

    def ufastqcplot(self,inline=True):
        """
        Return a mini-FastQC summary plot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        return ufastqcplot(self.fastqc.summary.path,inline=inline)

    def uscreenplot(self,inline=True):
        """
        Return a mini-FastQScreen summary plot

        Arguments:
          inline (bool): if True then return plot in format for
            inlining in HTML document
        """
        screen_files = list()
        for name in self.fastq_screen.names:
            screen_files.append(self.fastq_screen[name].txt)
        return uscreenplot(screen_files,inline=inline)

#######################################################################
# Functions
#######################################################################

def verify(project,qc_dir=None,qc_protocol=None):
    """
    Check the QC outputs are correct for a project

    Arguments:
      project (AnalysisProject): project to verify QC for
      qc_dir (str): path to the QC output dir; relative
        path will be treated as a subdirectory of the
        project being checked.
      qc_protocol (str): QC protocol to verify against
        (optional)

     Returns:
       Boolean: Returns True if all expected QC products
         are present, False if not.
    """
    logger.debug("verify: qc_dir (initial): %s" % qc_dir)
    if qc_dir is None:
        qc_dir = project.qc_dir
    else:
        if not os.path.isabs(qc_dir):
            qc_dir = os.path.join(project.dirn,
                                  qc_dir)
    logger.debug("verify: qc_dir (final)  : %s" % qc_dir)
    for dirn in (project.dirn,qc_dir):
        fastq_strand_conf = os.path.join(dirn,"fastq_strand.conf")
        if os.path.exists(fastq_strand_conf):
            break
        fastq_strand_conf = None
    logger.debug("verify: fastq_strand conf file : %s" %
                 fastq_strand_conf)
    cellranger_refdata = None
    qc_info_file = os.path.join(qc_dir,"qc.info")
    if os.path.exists(qc_info_file):
        qc_info = AnalysisProjectQCDirInfo(filen=qc_info_file)
        try:
            cellranger_refdata = qc_info['cellranger_refdata']
        except KeyError:
            pass
    logger.debug("verify: cellranger reference data : %s" %
                 cellranger_refdata)
    verified = True
    for f in expected_outputs(project,qc_dir,
                              fastq_strand_conf=fastq_strand_conf,
                              cellranger_refdata=cellranger_refdata,
                              qc_protocol=qc_protocol):
        if not os.path.exists(f):
            logging.debug("Missing: %s" % f)
            verified = False
    return verified

def report(projects,title=None,filename=None,qc_dir=None,
           report_attrs=None,summary_fields=None,
           relative_links=False,use_data_dir=False,
           make_zip=False):
    """
    Report the QC for a project

    Arguments:
      projects (list): AnalysisProject instances to report
        QC for
      title (str): optional, specify title for the report
        (defaults to '<PROJECT_NAME>: QC report')
      filename (str): optional, specify path and name for
        the output report file (defaults to
        '<PROJECT_NAME>.qc_report.html')
      qc_dir (str): path to the QC output dir
      report_attrs (list): optional, list of elements to
        report for each Fastq pair
      summary_fields (list): optional, list of fields to
        report for each sample in the summary table
        relative_links (boolean): optional, if set to True
        then use relative paths for links in the report
        (default is to use absolute paths)
      use_data_dir (boolean): if True then copy QC
        artefacts to a data directory parallel to the
        output report file
      make_zip (boolean): if True then also create a ZIP
        archive of the QC report and outputs (default is
        not to create the ZIP archive)

    Returns:
      String: filename of the output HTML report.
    """
    # Primary project
    project = projects[0]
    # Set output destination
    if filename is None:
        filename = "%s.qc_report.html" % project.name
    # Use data directory
    if use_data_dir:
        data_dir = "%s_data" % os.path.splitext(filename)[0]
    else:
        data_dir = None
    # Use relative paths for links
    if relative_links:
        relpath = os.path.dirname(filename)
    else:
        relpath = None
    # Initialise report
    report = QCReport(projects,
                      title=title,
                      qc_dir=qc_dir,
                      report_attrs=report_attrs,
                      summary_fields=summary_fields,
                      data_dir=data_dir,
                      relpath=relpath)
    # Styles
    report.add_css_rule(QC_REPORT_CSS_STYLES)
    # Write the report
    report.write(filename)
    # Make ZIP file
    if make_zip:
        # Name for ZIP file
        out_dir = os.path.dirname(filename)
        basename = os.path.splitext(os.path.basename(filename))[0]
        run_name = project.info.run
        zip_prefix = "%s.%s%s" % (basename,
                                  project.name,
                                  '.%s' % run_name if run_name else '')
        logging.debug("ZIP prefix: %s" % zip_prefix)
        zip_name = os.path.join(out_dir,"%s.zip" % zip_prefix)
        logging.debug("ZIP file: %s" % zip_name)
        # Directory with QC artefacts
        if data_dir:
            qc_dir = data_dir
        elif qc_dir is None:
            qc_dir = project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(project.dirn,
                                      qc_dir)
        logging.debug("report: qc_dir  = %s" % qc_dir)
        # Create ZIP archive
        try:
            tmp_zip = "%s.tmp" % zip_name
            zip_file = ZipArchive(tmp_zip,
                                  relpath=os.path.dirname(qc_dir),
                                  prefix=zip_prefix)
            # Add the report
            zip_file.add_file(filename)
            # Add the QC outputs
            logging.debug("Adding QC outputs to ZIP archive '%s'" % zip_name)
            for f in report.output_files:
                # Output files are absolute paths
                logging.debug("-- Adding %s" % f)
                if os.path.exists(f):
                    zip_file.add(f)
                else:
                    raise Exception("Unable to add missing file '%s'" % f)
            # Rename temporary ZIP to final location
            os.rename(tmp_zip,zip_name)
        except Exception as ex:
            # Remove temporary ZIP
            os.remove(tmp_zip)
            raise Exception("Failed to create ZIP archive '%s': %s" %
                            (zip_name,ex))
    # Return the output filename
    return filename

def pretty_print_reads(n):
    """
    Print the number of reads with commas at each thousand

    For example:

    >>> pretty_print_reads(10409789)
    10,409,789

    Arguments:
      n (int): number of reads

    Returns:
      String: representation with commas for every thousand.
    """
    n = str(int(n))[::-1]
    n0 = []
    while len(n) >= 3:
        n0.append(n[0:3])
        n = n[3:]
    if n: n0.append(n)
    return (','.join(n0))[::-1]

def sanitize_name(s,new_char='_'):
    """
    Replace 'unsafe' characters in HTML link targets

    Replaces 'unsafe' characters in a string used as a
    link target in an HTML document with underscore
    characters.

    Arguments:
      s (str): string to sanitize
      replace_with (str): string to replace 'unsafe'
        characters with (default: '_')
    """
    for c in ('#',':','/'):
        s = s.replace(c,new_char)
    return s
