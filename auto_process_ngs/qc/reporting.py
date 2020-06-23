#!/usr/bin/env python
#
#     reporting: report QC from analysis projects
#     Copyright (C) University of Manchester 2018-2019 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import logging
import time
import ast
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
        self._parent_dir = os.path.dirname(self._project.dirn)

    @property
    def name(self):
        return self._project.name

    @property
    def paired_end(self):
        return self._project.info.paired_end

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
        logger.debug("QCReporter.verify: qc_dir (initial): %s" % qc_dir)
        if qc_dir is None:
            qc_dir = self._project.qc_dir
        else:
            if not os.path.isabs(qc_dir):
                qc_dir = os.path.join(self._project.dirn,
                                      qc_dir)
        logger.debug("QCReporter.verify: qc_dir (final)  : %s" % qc_dir)
        for dirn in (self._project.dirn,qc_dir):
            fastq_strand_conf = os.path.join(dirn,"fastq_strand.conf")
            if os.path.exists(fastq_strand_conf):
                break
            fastq_strand_conf = None
        logger.debug("QCReporter.verify: fastq_strand conf file : %s" %
                     fastq_strand_conf)
        cellranger_refdata = None
        qc_info_file = os.path.join(qc_dir,"qc.info")
        if os.path.exists(qc_info_file):
            qc_info = AnalysisProjectQCDirInfo(filen=qc_info_file)
            try:
                cellranger_refdata = qc_info['cellranger_refdata']
            except KeyError:
                pass
        logger.debug("QCReporter.verify: cellranger reference data : %s" %
                     cellranger_refdata)
        verified = True
        for f in expected_outputs(self._project,qc_dir,
                                  fastq_strand_conf=fastq_strand_conf,
                                  cellranger_refdata=cellranger_refdata,
                                  qc_protocol=qc_protocol):
            if not os.path.exists(f):
                print("Missing: %s" % f)
                verified = False
        return verified

    def report(self,title=None,filename=None,qc_dir=None,
               qc_protocol=None,report_attrs=None,
               summary_fields=None,relative_links=False):
        """
        Report the QC for the project

        Arguments:
          title (str): optional, specify title for the report
            (defaults to '<PROJECT_NAME>: QC report')
          filename (str): optional, specify path and name for
            the output report file (defaults to
            '<PROJECT_NAME>.qc_report.html')
          qc_dir (str): path to the QC output dir
          qc_protocol (str): QC protocol to report against
            (optional)
          report_attrs (list): optional, list of elements to
            report for each Fastq pair
          summary_fields (list): optional, list of fields to
            report for each sample in the summary table
          relative_links (boolean): optional, if set to True
            then use relative paths for links in the report
            (default is to use absolute paths)

        Returns:
          String: filename of the output HTML report.
        """
        # Set title and output destination
        if title is None:
            title = "%s: QC report" % self.name
        if filename is None:
            filename = "%s.qc_report.html" % self.name
        # Use relative paths for links
        if relative_links:
            relpath = os.path.dirname(filename)
        else:
            relpath = None
        # Initialise report
        report = QCReport(self._project,
                          title=title,
                          qc_dir=qc_dir,
                          qc_protocol=qc_protocol,
                          report_attrs=report_attrs,
                          summary_fields=summary_fields,
                          relpath=relpath)
        # Styles
        report.add_css_rule(QC_REPORT_CSS_STYLES)
        # Write the report
        report.write(filename)
        # Return the output filename
        return filename

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
    Create a QC report document for a project

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
    def __init__(self,project,title=None,qc_dir=None,qc_protocol=None,
                 report_attrs=None,summary_fields=None,relpath=None):
        """
        Create a new QCReport instance

        Arguments:
          project (AnalysisProject): project to report QC for
          title (str): title for the report (defaults to
            "QC report: <PROJECT_NAME")
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
          qc_protocol (str): QC protocol to report against
            (optional)
          report_attrs (list): list of elements to report for
            each Fastq pair
          summary_fields (list): list of fields to report for
            each sample in the summary table
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
        logger.debug("QCReport: qc_dir (initial): %s" % qc_dir)
        # Status of report
        self.status = True
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
        logger.debug("QCReport: qc_dir (final): %s" % self.qc_dir)
        # How to read Fastq names
        self.fastq_attrs = self.project.fastq_attrs
        # Detect outputs
        self._detect_outputs()
        if self.samples:
            print("Samples found:")
            for sample in self.samples:
                print("\t- %s" % sample)
        else:
            logger.warning("No samples found")
        if self.fastqs:
            print("Fastqs referenced:")
            for fastq in self.fastqs:
                print("\t- %s" % fastq)
            print("Reads found:")
            for read in self.reads:
                print("\t- %s" % read)
        else:
            logger.warning("No Fastqs referenced")
        if self.outputs:
            print("Available QC outputs:")
            for output in self.outputs:
                print("\t- %s" % output)
        else:
            logger.warning("No QC outputs found")
        if self.software:
            print("Software versions:")
            for package in self.software:
                print("\t- %s: %s" % (package,
                                      ','.join(self.software[package])))
        # Load metadata from parent run
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
        # Set up title
        if title is None:
            title = "QC report: %s" % self.project.name
        # Initialise superclass
        Document.__init__(self,title)
        # Relative paths
        if relpath is not None:
            relpath = os.path.normpath(os.path.abspath(relpath))
        self.relpath = relpath
        # Field descriptions for summary table
        self.field_descriptions = {
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
        # Fields to report in summary table
        if not summary_fields:
            if len(self.reads) > 1:
                summary_fields = ['sample',
                                  'fastqs',
                                  'reads',
                                  'read_lengths',]
            else:
                summary_fields = ['sample',
                                  'fastq',
                                  'reads',
                                  'read_lengths',]
            if 'strandedness' in self.outputs:
                summary_fields.append('strandedness')
            for read in self.reads:
                if ('fastqc_%s' % read) in self.outputs:
                    summary_fields.append('fastqc_%s' % read)
                    summary_fields.append('boxplot_%s' % read)
                if ('screens_%s' % read) in self.outputs:
                    summary_fields.append('screens_%s' % read)
            if 'cellranger_count' in self.outputs:
                summary_fields.append('cellranger_count')
        self.summary_fields = summary_fields
        # Attributes to report for each sample
        if report_attrs is None:
            report_attrs = ['fastqc',
                            'fastq_screen',]
            if 'strandedness' in self.outputs:
                report_attrs.append('strandedness')
        self.report_attrs = report_attrs
        # Initialise tables
        self.metadata_table = self._init_metadata_table()
        self.software_table = self._init_software_table()
        self.summary_table = self._init_summary_table()
        # Initialise report sections
        self.preamble = self._init_preamble_section()
        self.warnings = self._init_warnings_section()
        self.summary = self._init_summary_section()
        # Build the report
        print("Building the report...")
        self.report_metadata()
        self.report_software()
        for sample in self.samples:
            self.report_sample(sample)
        # Report the status
        self.report_status()

    def _init_metadata_table(self):
        """
        Internal: set up a table for project metadata

        Associated CSS class is 'metadata'
        """
        metadata_tbl = Table(('item','value',))
        metadata_tbl.no_header()
        metadata_tbl.add_css_classes('metadata')
        return metadata_tbl

    def _init_software_table(self):
        """
        Internal: set up a table for software information

        Associated CSS class is 'metadata'
        """
        software_tbl = Table(('program','version',))
        software_tbl.no_header()
        software_tbl.add_css_classes('metadata')
        return software_tbl

    def _init_summary_table(self):
        """
        Internal: set up a table for summarising samples

        Associated CSS classes are 'summary' and 'fastq_summary'
        """
        summary_tbl = Table(self.summary_fields,
                            **self.field_descriptions)
        summary_tbl.add_css_classes('summary','fastq_summary')
        return summary_tbl

    def _init_preamble_section(self):
        """
        Internal: set up a "preamble" section
        """
        preamble = self.add_section()
        preamble.add("Report generated by auto_process %s on %s" %
                     (get_version(),time.asctime()))
        return preamble

    def _init_summary_section(self):
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
        comments = info.add_subsection("Comments",
                                       css_classes=("info",))
        comments_list = List()
        try:
            if self.project.info.comments:
                for comment in self.project.info.comments.split(';'):
                    comments_list.add_item(comment.strip())
            else:
                # Drop out with exception
                raise AttributeError
        except AttributeError:
            comments_list.add_item("N/A")
        comments.add(comments_list)
        # Add an empty section to clear HTML floats
        clear = summary.add_subsection(css_classes=("clear",))
        # Add the summary table
        summary.add("%d samples | %d fastqs" % (len(self.project.samples),
                                                len(self.project.fastqs)))
        summary.add(self.summary_table)
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

    def _detect_outputs(self):
        """
        Internal: determine which QC outputs are present
        """
        outputs = set()
        software = {}
        print("Scanning contents of %s" % self.qc_dir)
        files = os.listdir(self.qc_dir)
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
                versions.add(Fastqscreen(
                    os.path.join(self.qc_dir,screen)).version)
            if versions:
                software['fastq_screen'] = sorted(list(versions))
        # Look for fastqc outputs
        fastqcs = list(filter(lambda f:
                              f.endswith("_fastqc.html"),
                              files))
        logger.debug("Fastqc: %s" % fastqcs)
        print("\t- %d fastqc files" % len(fastqcs))
        if fastqcs:
            versions = set()
            # Pull out the Fastq names from the Fastqc files
            for fastqc in fastqcs:
                fastqc = os.path.splitext(fastqc)[0]
                f = os.path.basename(fastqc)[:-len("_fastqc")]
                fastq_names.add(f)
                fq = self.fastq_attrs(f)
                outputs.add("fastqc_%s%s" %
                            (('i' if fq.is_index_read else 'r'),
                             (fq.read_number
                              if fq.read_number is not None else '1')))
                versions.add(Fastqc(
                    os.path.join(self.qc_dir,fastqc)).version)
            if versions:
                software['fastqc'] = sorted(list(versions))
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
                versions.add(Fastqstrand(
                    os.path.join(self.qc_dir,f)).version)
            if versions:
                software['fastq_strand'] = sorted(list(versions))
        # Look for ICELL8 outputs
        print("Checking for ICELL8 reports in %s/stats" %
              self.project.dirn)
        icell8_stats_xlsx = os.path.join(self.project.dirn,
                                         "stats",
                                         "icell8_stats.xlsx")
        if os.path.exists(icell8_stats_xlsx):
            outputs.add("icell8_stats")
        icell8_report_html = os.path.join(self.project.dirn,
                                          "icell8_processing.html")
        if os.path.exists(icell8_report_html):
            outputs.add("icell8_report")
        # Look for cellranger_count outputs
        cellranger_count_dir = os.path.join(self.qc_dir,
                                            "cellranger_count")
        cellranger_samples = []
        if os.path.isdir(cellranger_count_dir):
            for d in filter(
                    lambda f:
                    os.path.isdir(os.path.join(cellranger_count_dir,f)),
                    os.listdir(cellranger_count_dir)):
                cellranger_samples.append(d)
            if cellranger_samples:
                outputs.add("cellranger_count")
        # Look for MultiQC report
        print("Checking for MultiQC report in %s" % self.project.dirn)
        multiqc_report = os.path.join(self.project.dirn,
                                      "multi%s_report.html"
                                      % os.path.basename(self.qc_dir))
        if os.path.isfile(multiqc_report):
            outputs.add("multiqc")
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
        return self.outputs

    def report_metadata(self):
        """
        Report the project metadata

        Adds entries for the project metadata to the "metadata"
        table in the report
        """
        metadata_items = ['run',
                          'run_id',
                          'user',
                          'PI',
                          'library_type',
                          'single_cell_platform',
                          'number_of_cells',
                          'organism',
                          'protocol',]
        if 'cellranger_count' in self.outputs:
            metadata_items.append('cellranger_reference')
        if 'multiqc' in self.outputs:
            metadata_items.append('multiqc')
        if 'icell8_stats' in self.outputs:
            metadata_items.append('icell8_stats')
        if 'icell8_report' in self.outputs:
            metadata_items.append('icell8_report')
        metadata_titles = {
            'run_id': 'Run ID',
            'run': 'Run name',
            'user': 'User',
            'PI': 'PI',
            'library_type': 'Library type',
            'single_cell_platform': 'Single cell preparation platform',
            'number_of_cells': 'Number of cells',
            'organism': 'Organism',
            'protocol': 'QC protocol',
            'cellranger_reference': 'Cellranger reference dataset',
            'multiqc': 'MultiQC report',
            'icell8_stats': 'ICELL8 statistics',
            'icell8_report': 'ICELL8 processing report',
        }
        for item in metadata_items:
            # Try to acquire the value from QC metadata
            try:
                value = self.project.qc_info(self.qc_dir)[item]
            except KeyError:
                # Fall back to project metadata
                try:
                    if self.project.info[item]:
                        value = self.project.info[item]
                    else:
                        # No value set, skip this item
                        continue
                except KeyError:
                    # Additional non-metadata items, or items
                    # requiring additional processing
                    if item == 'run_id':
                        try:
                            value = run_reference_id(
                                self.project.info['run'],
                                platform=self.project.info['platform'],
                                facility_run_number=
                                self.run_metadata['run_number'])
                        except (AttributeError,TypeError) as ex:
                            logger.warning("Run reference ID can't be "
                                           "determined: %s (ignored)" % ex)
                            continue
                    elif item == 'cellranger_reference':
                        path = self.project.qc_info(self.qc_dir).\
                               cellranger_refdata
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
                                         % os.path.basename(self.qc_dir)
                        value = Link(multiqc_report)
                    elif item == 'icell8_stats':
                        value = Link("icell8_stats.xlsx",
                                     os.path.join("stats",
                                                  "icell8_stats.xlsx"))
                    elif item == 'icell8_report':
                        value = Link("icell8_processing.html")
                    else:
                        raise Exception("Unrecognised item to report: '%s'"
                                        % item)
            # Add to the metadata table
            self.metadata_table.add_row(
                item=metadata_titles[item],
                value=value)

    def report_software(self):
        """
        Report the software versions used in the processing & QC

        Adds entries for the software versions to the "software"
        table in the report
        """
        software_packages = ['bcl2fastq',
                             'cellranger',
                             'fastqc',
                             'fastq_screen',
                             'fastq_strand',]
        software_names = {
            'bcl2fastq': 'Bcl2fastq',
            'cellranger': 'Cellranger',
            'fastqc': 'FastQC',
            'fastq_screen': 'FastqScreen',
            'fastq_strand': 'FastqStrand',
        }
        for pkg in software_packages:
            # Acquire the value
            try:
                if self.software[pkg]:
                    value = ','.join(self.software[pkg])
                else:
                    # No value set, skip this item
                    continue
            except KeyError:
                if pkg == 'bcl2fastq':
                    try:
                        bcl2fastq = ast.literal_eval(
                            self.run_metadata.bcl2fastq_software)
                        value = bcl2fastq[2]
                    except ValueError:
                        continue
                elif pkg == 'cellranger':
                    try:
                        cellranger = ast.literal_eval(
                            self.run_metadata.cellranger_software)
                        software_names['cellranger'] = cellranger[1].title()
                        value = cellranger[2]
                    except ValueError:
                        continue
                elif pkg not in software_names:
                    # Unrecognised package name
                    raise Exception("Unrecognised software package "
                                    "item to report: '%s'" % pkg)
                else:
                    # No value set, skip this item
                    continue
            # Add to the metadata table
            self.software_table.add_row(
                program=software_names[pkg],
                version=value)

    def report_sample(self,sample):
        """
        Report the QC for a sample

        Reports the QC for the sample and Fastqs to
        the summary table and appends a section with
        detailed reports to the document.

        Arguments:
          sample (str): name of sample to report
        """
        sample_report = self.add_section(
            "Sample: %s" % sample,
            name="sample_%s" % sample,
            css_classes=('sample',))
        # Get Fastq groups
        fastqs = sorted(list(
            filter(lambda fq:
                   self.fastq_attrs(fq).sample_name == sample,
                   self.fastqs)))
        fastq_groups = group_fastqs_by_name(fastqs,self.fastq_attrs)
        # Number of fastqs
        if len(self.reads) > 1:
            sample_report.add("%d fastq R1/R2 pairs" %
                              len(fastq_groups))
        else:
            sample_report.add("%d fastqs" %
                              len(fastq_groups))
        # Keep track of the first line in the summary
        # table, as per-sample metrics (and name)
        # should only be reported on the first line
        first_line = True
        # Report each Fastq group
        for fastqs in fastq_groups:
            # Report Fastq pair
            fastq_group = QCReportFastqGroup(fastqs,
                                             qc_dir=self.qc_dir,
                                             fastq_attrs=self.fastq_attrs)
            fastq_group.report(sample_report,
                               attrs=self.report_attrs,
                               relpath=self.relpath)
            # Add line in summary table
            if sample is not None:
                if first_line:
                    # Only display sample name on first line
                    idx = self.summary_table.add_row(
                        sample=Link(sample,
                                    sample_report))
                else:
                    # Don't display sample name for subsequent
                    # lines in summary table
                    idx = self.summary_table.add_row(sample="&nbsp;")
            else:
                idx = self.summary_table.add_row(sample="&nbsp;")
            status = fastq_group.update_summary_table(
                self.summary_table,idx=idx,
                fields=self.summary_fields,
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
    def __init__(self,fastqs,qc_dir,fastq_attrs=AnalysisFastq):
        """
        Create a new QCReportFastqGroup

        Arguments:
          fastqs (list): list of paths to Fastqs in the
            group
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
          fastq_attrs (BaseFastqAttrs): class for extracting
            data from Fastq names
        """
        self.qc_dir = qc_dir
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
        Add report for Fastq pair to a document section

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
                        print("Dropping %s" % field)
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
                logger.error("Exception setting '%s' in summary table "
                              "for Fastq group { %s }: %s (ignored)"
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
    def __init__(self,fastq,qc_dir,fastq_attrs=AnalysisFastq):
        """
        Create a new QCReportFastq instance

        Arguments:
          fastq (str): path to Fastq file
          qc_dir (str): path to QC directory
          fastq_attrs (BaseFastqAttrs): class for extracting
            data from Fastq names
        """
        # Source data
        self.name = os.path.basename(fastq)
        self.path = os.path.abspath(fastq)
        self.safe_name = strip_ngs_extensions(self.name)
        self.fastq_attrs = fastq_attrs
        # Sample name
        self.sample_name = fastq_attrs(self.name).sample_name
        # FastQC
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
