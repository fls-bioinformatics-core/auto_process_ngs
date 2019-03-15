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
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import cmp_sample_names
from bcftbx.TabFile import TabFile
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from ..analysis import run_reference_id
from ..docwriter import Document
from ..docwriter import Section
from ..docwriter import Table
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import Target
from ..docwriter import List
from ..metadata import AnalysisDirMetadata
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .fastq_strand import Fastqstrand
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
.clear { clear: both; }
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
table.summary td { text-align: right;
                   padding: 2px 5px;
                   border-bottom: solid 1px lightgray; }
table.fastq_summary tr td:first-child {
          background-color: grey;
          color: white;
          font-weight: bold; }
table.fastq_summary tr td:first-child a {
          color: white;
          font-weight: bold; }
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
    samples: list of QCSample instances

    Provides the follow methods:

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
        self._samples = []
        self._parent_dir = os.path.dirname(self._project.dirn)
        for sample in self._project.samples:
            self._samples.append(QCSample(sample))
        self._samples = sorted(self._samples,
                               cmp=lambda x,y: cmp_sample_names(x.name,y.name))
        logger.debug("Found %d samples" % len(self._samples))

    @property
    def name(self):
        return self._project.name

    @property
    def paired_end(self):
        return self._project.info.paired_end

    @property
    def samples(self):
        return self._samples

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
        fastq_strand_conf = os.path.join(self._project.dirn,
                                         "fastq_strand.conf")
        logger.debug("QCReporter.verify: qc_dir (final)  : %s" % qc_dir)
        verified = True
        for f in expected_outputs(self._project,qc_dir,
                                  fastq_strand_conf,
                                  qc_protocol=qc_protocol):
            if not os.path.exists(f):
                print "Missing: %s" % f
                verified = False
                ##break
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

class QCSample(object):
    """
    Class describing QC results for an AnalysisSample

    Provides the follow properties:

    name: sample name
    fastq_pairs: list of FastqSet instances

    Provides the follow methods:

    verify: checks that QC outputs are present
    """
    def __init__(self,sample):
        """
        Initialise a new QCSample instance

        Arguments:
           sample (AnalysisSample): sample instance
        """
        self._sample = sample
        self._fastq_pairs = get_fastq_pairs(sample)

    @property
    def name(self):
        return self._sample.name

    @property
    def fastq_pairs(self):
        return self._fastq_pairs

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
        return filter(lambda fq: fq is not None,
                      self._fastqs)

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
    - fastqc_r1: FastQC mini-plot for R1
    - boxplot_r1: FastQC per-base-quality mini-boxplot' for R1
    - screens_r1: FastQScreen mini-plots for R1
    - fastqc_r2: FastQC mini-plot for R2
    - boxplot_r2: FastQC per-base-quality mini-boxplot' for R2
    - screens_r2: FastQScreen mini-plots for R2
    - strandedness: 'forward', 'reverse' or 'unstranded' for pair

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
        # Detect outputs
        self._detect_outputs()
        if self.outputs:
            print "Available QC outputs:"
            for output in self.outputs:
                print "\t- %s" % output
        else:
            logger.warning("No QC outputs found")
        if self.software:
            print "Software versions:"
            for package in self.software:
                print "\t- %s: %s" % (package,
                                      ','.join(self.software[package]))
        # Load metadata from parent run
        self.run_metadata = AnalysisDirMetadata()
        run_metadata_file = os.path.join(
            os.path.dirname(os.path.abspath(self.project.dirn)),
            "metadata.info")
        if os.path.exists(run_metadata_file):
            print "Loading run metadata from %s" % run_metadata_file
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
        self.field_descriptions = { 'sample': 'Sample',
                                    'fastq' : 'Fastq',
                                    'fastqs': 'Fastqs (R1/R2)',
                                    'reads': '#reads',
                                    'read_lengths': 'Length',
                                    'fastqc_r1': 'FastQC[R1]',
                                    'boxplot_r1': 'Boxplot[R1]',
                                    'screens_r1': 'Screens[R1]',
                                    'fastqc_r2': 'FastQC[R2]',
                                    'boxplot_r2': 'Boxplot[R2]',
                                    'screens_r2': 'Screens[R2]',
                                    'strandedness': 'Strand', }
        # Fields to report in summary table
        if not summary_fields:
            if self.project.info.paired_end:
                reads = ('r1','r2')
                summary_fields = ['sample',
                                  'fastqs',
                                  'reads',
                                  'read_lengths',]
            else:
                reads = ('r1',)
                summary_fields = ['sample',
                                  'fastq',
                                  'reads',
                                  'read_lengths',]
            if 'strandedness' in self.outputs:
                summary_fields.append('strandedness')
            for read in reads:
                if ('fastqc_%s' % read) in self.outputs:
                    summary_fields.append('fastqc_%s' % read)
                    summary_fields.append('boxplot_%s' % read)
                if ('screens_%s' % read) in self.outputs:
                    summary_fields.append('screens_%s' % read)
        self.summary_fields = summary_fields
        # Attributes to report for each sample
        if report_attrs is None:
            report_attrs = ['fastqc',
                            'fastq_screen',
                            'program_versions']
            if 'strandedness' in self.outputs:
                report_attrs.append('strandedness')
        self.report_attrs = report_attrs
        # Initialise tables
        self.metadata_table = self._init_metadata_table()
        self.software_table = self._init_software_table()
        self.summary_table = self._init_summary_table()
        # Initialise report sections
        self.preamble = self._init_preamble_section()
        self.summary = self._init_summary_section()
        # Build the report
        print "Building the report..."
        self.report_metadata()
        self.report_software()
        for sample in self.project.samples:
            try:
                self.report_sample(sample)
            except Exception as ex:
                logger.error("Exception for sample '%s': %s"
                             % (sample.name,ex))
                raise ex

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

    def _detect_outputs(self):
        """
        Internal: determine which QC outputs are present
        """
        outputs = []
        software = {}
        print "Scanning contents of %s" % self.qc_dir
        files = os.listdir(self.qc_dir)
        print "\t- %d objects found" % len(files)
        fastqs = [strip_ngs_extensions(os.path.basename(fq))
                  for fq in self.project.fastqs]
        fastqs_r1 = filter(lambda f:
                           self.project.fastq_attrs(f).read_number == 1,
                           fastqs)
        fastqs_r2 = filter(lambda f:
                           self.project.fastq_attrs(f).read_number == 2,
                           fastqs)
        logger.debug("files: %s" % files)
        logger.debug("fastqs: %s" % fastqs)
        logger.debug("fastqs R1: %s" % fastqs_r1)
        logger.debug("fastqs R2: %s" % fastqs_r2)
        # Look for screen files
        screens = filter(lambda f:
                         f.endswith("_screen.txt") or
                         f.endswith("_screen.png"),
                         files)
        logger.debug("Screens: %s" % screens)
        print "\t- %d fastq_screen files" % len(screens)
        if screens:
            for fq in fastqs_r1:
                if filter(lambda s: s.startswith(fq),screens):
                    outputs.append("screens_r1")
                    break
            for fq in fastqs_r2:
                if filter(lambda s: s.startswith(fq),screens):
                    outputs.append("screens_r2")
                    break
            versions = set()
            for screen in filter(lambda s:
                                 s.endswith("_screen.txt"),
                                 screens):
                versions.add(Fastqscreen(
                    os.path.join(self.qc_dir,screen)).version)
            if versions:
                software['fastq_screen'] = sorted(list(versions))
        # Look for fastqc outputs
        fastqcs = filter(lambda f: f.endswith("_fastqc.html"),files)
        logger.debug("Fastqc: %s" % fastqcs)
        print "\t- %d fastqc files" % len(fastqcs)
        if fastqcs:
            for fq in fastqs_r1:
                if filter(lambda f: f.startswith(fq),fastqcs):
                    outputs.append("fastqc_r1")
                    break
            for fq in fastqs_r2:
                if filter(lambda f: f.startswith(fq),fastqcs):
                    outputs.append("fastqc_r2")
                    break
            versions = set()
            for f in fastqcs:
                d = os.path.splitext(f)[0]
                versions.add(Fastqc(
                    os.path.join(self.qc_dir,d)).version)
            if versions:
                software['fastqc'] = sorted(list(versions))
        # Look for fastq_strand outputs
        fastq_strand = filter(lambda f: f.endswith("_fastq_strand.txt"),files)
        logger.debug("fastq_strand: %s" % fastq_strand)
        print "\t- %d fastq_strand files" % len(fastq_strand)
        if fastq_strand:
            outputs.append("strandedness")
            versions = set()
            for f in fastq_strand:
                versions.add(Fastqstrand(
                    os.path.join(self.qc_dir,f)).version)
            if versions:
                software['fastq_strand'] = sorted(list(versions))
        # Look for MultiQC report
        print "Checking for MultiQC report in %s" % self.project.dirn
        multiqc_report = os.path.join(self.project.dirn,
                                      "multi%s_report.html"
                                      % os.path.basename(self.qc_dir))
        if os.path.isfile(multiqc_report):
            outputs.append("multiqc")
        self.outputs = outputs
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
                          'qc_protocol',]
        if 'multiqc' in self.outputs:
            metadata_items.append('multiqc')
        metadata_titles = {
            'run_id': 'Run ID',
            'run': 'Run name',
            'user': 'User',
            'PI': 'PI',
            'library_type': 'Library type',
            'single_cell_platform': 'Single cell preparation platform',
            'number_of_cells': 'Number of cells',
            'organism': 'Organism',
            'qc_protocol': 'QC protocol',
            'multiqc': 'MultiQC report',
        }
        for item in metadata_items:
            # Acquire the value
            try:
                if self.project.info[item]:
                    value = self.project.info[item]
                else:
                    # No value set, skip this item
                    continue
            except KeyError:
                if item == 'run_id':
                    try:
                        value = run_reference_id(
                            self.project.info['run'],
                            platform=self.project.info['platform'],
                            facility_run_number=
                            self.run_metadata['run_number'])
                    except AttributeError as ex:
                        logger.warning("Run reference ID can't be "
                                       "determined: %s (ignored)" % ex)
                        continue
                elif item == 'qc_protocol':
                    value = self.project.qc_info(self.qc_dir).protocol
                elif item == 'multiqc':
                    multiqc_report = "multi%s_report.html" \
                                     % os.path.basename(self.qc_dir)
                    value = Link(multiqc_report)
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
          sample (AnalysisSample): sample to report
        """
        sample = QCSample(sample)
        # Create a new section for the sample
        sample_report = self.add_section(
            "Sample: %s" % sample.name,
            name="sample_%s" % sample.name,
            css_classes=('sample',))
        # Number of fastqs
        if self.project.info.paired_end:
            sample_report.add("%d fastq R1/R2 pairs" %
                              len(sample.fastq_pairs))
        else:
            sample_report.add("%d fastqs" %
                              len(sample.fastq_pairs))
        # Report each Fastq/Fastq pair
        sample_name = sample.name
        for fq_pair in sample.fastq_pairs:
            # Report Fastq pair
            fastq_pair = QCReportFastqPair(fq_pair.r1,
                                           fq_pair.r2,
                                           self.qc_dir)
            fastq_pair.report(sample_report,
                              attrs=self.report_attrs,
                              relpath=self.relpath)
            # Add line in summary table
            if sample_name is not None:
                idx = self.summary_table.add_row(sample=Link(sample_name,
                                                             sample_report))
            else:
                idx = self.summary_table.add_row(sample="&nbsp;")
            fastq_pair.update_summary_table(self.summary_table,idx=idx,
                                            fields=self.summary_fields,
                                            relpath=self.relpath)

class QCReportFastqPair(object):
    """
    Utility class for reporting the QC for a Fastq pair

    Provides the following properties:

    r1: QCReportFastq instance for R1 Fastq
    r2: QCReportFastq instance for R2 Fastq

    Provides the following methods:

    report: 
    """
    def __init__(self,fastqr1,fastqr2,qc_dir):
        """
        Create a new QCReportFastqPair

        Arguments:
          fastqr1 (str): R1 Fastq file
          fastqr2 (str): R2 Fastq file (None if 'pair' is
            single-ended)
          qc_dir (str): path to the QC output dir; relative
            path will be treated as a subdirectory of the
            project
        """
        self.fastqr1 = fastqr1
        self.fastqr2 = fastqr2
        self.qc_dir = qc_dir

    @property
    def paired_end(self):
        """
        True if pair consists of R1/R2 files
        """
        return (self.fastqr2 is not None)

    @property
    def r1(self):
        """
        QCReportFastq instance for R1 Fastq
        """
        return QCReportFastq(self.fastqr1,self.qc_dir)

    @property
    def r2(self):
        """
        QCReportFastq instance for R2 Fastq (None if not paired end)
        """
        if self.fastqr2 is not None:
            return QCReportFastq(self.fastqr2,self.qc_dir)
        return None

    @property
    def fastq_strand_txt(self):
        """
        Locate output from fastq_strand (None if not found)
        """
        fastq_strand_txt = None
        for fq in (self.fastqr1,self.fastqr2):
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
            name="strandedness_%s" % self.r1.safe_name)
        strandedness_report.add_css_classes("strandedness")
        # Locate strandedness
        txt = self.fastq_strand_txt
        # No file found
        if txt is None:
            strandedness_report.add("!!!No strandedness data available!!!")
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
        for fq in (self.r1,self.r2):
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
                    fq.report_program_versions(fq_report)
        # Sections for pairwise data
        if "strandedness" in attrs:
            # Strandedness
            self.report_strandedness(fastqs_report)
        # Add an empty section to clear HTML floats
        clear = fastqs_report.add_subsection(css_classes=("clear",))

    def update_summary_table(self,summary_table,idx=None,fields=None,
                             relpath=None):
        """
        Add a line to a summary table reporting a Fastq pair

        Creates a new line in 'summary_table' (or updates an
        existing line) for the Fastq pair, adding content for
        each specied field.

        The following fields can be reported for each Fastq
        pair:

        - fastqs (if paired-end)
        - fastq (if single-end)
        - reads
        - read_lengths
        - boxplot_r1
        - boxplot_r2
        - fastqc_r1
        - fastqc_r2
        - screens_r1
        - screens_r2
        - strandedness

        Arguments:
          summary_table (Table): table to add the summary to
          idx (int): if supplied then indicates which existing
            table row to update (if None then a new row is
            appended)
          fields (list): list of custom fields to report
          relpath (str): if set then make link paths
            relative to 'relpath'
        """
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
        # Add row to summary table
        if idx is None:
            idx = summary_table.add_row()
        # Populate with data
        for field in fields:
            try:
                if field == "sample":
                    logger.debug("'sample' ignored")
                    continue
                elif field == "fastq":
                    value = Link(self.r1.name,"#%s" % self.r1.safe_name)
                elif field == "fastqs":
                    value = "%s<br />%s" % \
                            (Link(self.r1.name,
                                  "#%s" % self.r1.safe_name),
                             Link(self.r2.name,
                                  "#%s" % self.r2.safe_name))
                elif field == "reads":
                    value = pretty_print_reads(
                        self.r1.fastqc.data.basic_statistics(
                            'Total Sequences'))
                elif field == "read_lengths":
                    read_lengths_r1 = Link(
                        self.r1.fastqc.data.basic_statistics(
                            'Sequence length'),
                        self.r1.fastqc.summary.link_to_module(
                            'Sequence Length Distribution',
                            relpath=relpath))
                    if self.paired_end:
                        read_lengths_r2 = Link(
                            self.r2.fastqc.data.basic_statistics(
                                'Sequence length'),
                            self.r2.fastqc.summary.link_to_module(
                                'Sequence Length Distribution',
                                relpath=relpath))
                        read_lengths = "%s<br />%s" % (read_lengths_r1,
                                                       read_lengths_r2)
                    else:
                        read_lengths = "%s" % read_lengths_r1
                    value = read_lengths
                elif field == "boxplot_r1":
                    value = Img(self.r1.uboxplot(),
                                href="#boxplot_%s" %
                                self.r1.safe_name)
                elif field == "boxplot_r2":
                    value = Img(self.r2.uboxplot(),
                                href="#boxplot_%s" %
                                self.r2.safe_name)
                elif field == "fastqc_r1":
                    value = Img(self.r1.ufastqcplot(),
                                href="#fastqc_%s" %
                                self.r1.safe_name,
                                title=self.r1.fastqc_summary())
                elif field == "fastqc_r2":
                    value = Img(self.r2.ufastqcplot(),
                                href="#fastqc_%s" %
                                self.r2.safe_name,
                                title=self.r2.fastqc_summary())
                elif field == "screens_r1":
                    value = Img(self.r1.uscreenplot(),
                                href="#fastq_screens_%s" %
                                self.r1.safe_name)
                elif field == "screens_r2":
                    value = Img(self.r2.uscreenplot(),
                                href="#fastq_screens_%s" %
                                self.r2.safe_name)
                elif field == "strandedness":
                    value = Img(self.ustrandplot(),
                                href="#strandedness_%s" %
                                self.r1.safe_name,
                                title=self.strandedness())
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
                if self.paired_end:
                    fastqs = "pair %s/%s" % (self.r1.name,
                                             self.r2.name)
                else:
                    fastqs = "%s" % self.r1.name
                logger.error("Exception setting '%s' in summary table "
                              "for Fastq %s: %s" % (field,
                                                    fastqs,
                                                    ex))
                # Put error value into the table
                summary_table.set_value(idx,field,"<b>ERROR</b>")

class QCReportFastq(object):
    """
    Provides interface to QC outputs for Fastq file

    Provides the following attributes:

    name: basename of the Fastq
    path: path to the Fastq
    safe_name: name suitable for use in HTML links etc
    fastqc: Fastqc instance
    fastq_screen.names: list of FastQScreen names
    fastq_screen.SCREEN.description: description of SCREEN
    fastq_screen.SCREEN.png: associated PNG file for SCREEN
    fastq_screen.SCREEN.txt: associated TXT file for SCREEN
    fastq_screen.SCREEN.version: associated version for SCREEN
    program_versions.NAME: version of package NAME

    Provides the following methods:

    report_fastqc
    report_fastq_screens
    report_program_versions
    uboxplot
    ufastqcplot
    uscreenplot
    """
    def __init__(self,fastq,qc_dir):
        """
        Create a new QCReportFastq instance

        Arguments:
          fastq (str): path to Fastq file
          qc_dir (str): path to QC directory
        """
        # Source data
        self.name = os.path.basename(fastq)
        self.path = os.path.abspath(fastq)
        self.safe_name = strip_ngs_extensions(self.name)
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
            fastqc_version = '?'
        self.program_versions['fastqc'] = fastqc_version
        fastq_screen_versions = list(
            set([self.fastq_screen[s].version
                 for s in self.fastq_screen.names]))
        if fastq_screen_versions:
            fastq_screen_versions = ','.join(sorted(fastq_screen_versions))
        else:
            fastq_screen_versions = '?'
        self.program_versions['fastq_screen'] = fastq_screen_versions

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
            fastqc_report.add("!!!No FastQC data available!!!")
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
            screens_report.add("No screens found")
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

def get_fastq_pairs(sample):
    """
    Return pairs of Fastqs for an AnalysisSample instance

    Arguments:
       sample (AnalysisSample): sample to get Fastq pairs for

    Returns:
       list: list of FastqSet instances, sorted by R1 names
    """
    pairs = []
    fastqs_r1 = sample.fastq_subset(read_number=1)
    fastqs_r2 = sample.fastq_subset(read_number=2)
    for fqr1 in fastqs_r1:
        # Split up R1 name
        logger.debug("fqr1 %s" % os.path.basename(fqr1))
        dir_path = os.path.dirname(fqr1)
        # Generate equivalent R2 file
        fqr2 = sample.fastq_attrs(fqr1)
        fqr2.read_number = 2
        fqr2 = os.path.join(dir_path,"%s%s" % (fqr2,fqr2.extension))
        logger.debug("fqr2 %s" % os.path.basename(fqr2))
        if fqr2 in fastqs_r2:
            pairs.append(FastqSet(fqr1,fqr2))
        else:
            pairs.append(FastqSet(fqr1))
    pairs = sorted(pairs,cmp=lambda x,y: cmp(x.r1,y.r1))
    return pairs

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
