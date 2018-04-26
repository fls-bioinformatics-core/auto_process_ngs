#!/usr/bin/env python
#
#     reporting: report QC from analysis projects
#     Copyright (C) University of Manchester 2018 Peter Briggs
#

#######################################################################
# Imports
#######################################################################

import sys
import os
import logging
import time
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import cmp_sample_names
from bcftbx.TabFile import TabFile
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import AttributeDictionary
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from ..applications import Command
from ..docwriter import Document
from ..docwriter import Section
from ..docwriter import Table
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import Target
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .illumina_qc import IlluminaQC
from .illumina_qc import fastqc_output
from .illumina_qc import fastq_screen_output
from .plots import uscreenplot
from .plots import ufastqcplot
from .plots import uboxplot
from .plots import encode_png
from .. import get_version

# Data
from .illumina_qc import FASTQ_SCREENS

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class QCReporter(object):
    """
    Class describing QC results for an AnalysisProject

    """
    def __init__(self,project):
        """
        Initialise a new QCReporter instance

        Arguments:
           project (AnalysisProject): project to handle the QC for

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

    def verify(self,qc_dir=None):
        """
        Check that the QC outputs are correct

        Returns True if the QC appears to have run successfully,
        False if not.

        Arguments:
          qc_dir (str): path to the QC output dir

        """
        if qc_dir is None:
            qc_dir = self._project.qc_dir
        verified = True
        for sample in self._samples:
            if not sample.verify(qc_dir):
                verified = False
        return verified

    def report(self,title=None,filename=None,qc_dir=None,
               relative_links=False):
        """
        Report the QC for the project

        Arguments:
          title (str): optional, specify title for the report
            (defaults to '<PROJECT_NAME>: QC report')
          filename (str): optional, specify path and name for
            the output report file (defaults to
            '<PROJECT_NAME>.qc_report.html')
          qc_dir (str): path to the QC output dir
          relative_links (boolean): optional, if set to True
            then use relative paths for links in the report
            (default is to use absolute paths)

        Returns:
          String: filename of the HTML report.

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
        report = QCReport(self._project,title=title)
        # Styles
        report.add_css_rule("h1 { background-color: #42AEC2;\n"
                            "     color: white;\n"
                            "     padding: 5px 10px; }")
        report.add_css_rule("h2 { background-color: #8CC63F;\n"
                            "     color: white;\n"
                            "     display: inline-block;\n"
                            "     padding: 5px 15px;\n"
                            "     margin: 0;\n"
                            "     border-top-left-radius: 20px;\n"
                            "     border-bottom-right-radius: 20px; }")
        report.add_css_rule("h3, h4 { background-color: grey;\n"
                            "     color: white;\n"
                            "     display: block;\n"
                            "     padding: 5px 15px;\n"
                            "     margin: 0;\n"
                            "     border-top-left-radius: 20px;\n"
                            "     border-bottom-right-radius: 20px; }")
        report.add_css_rule(".sample { margin: 10 10;\n"
                            "          border: solid 2px #8CC63F;\n"
                            "          padding: 0;\n"
                            "          background-color: #ffe;\n"
                            "          border-top-left-radius: 25px;\n"
                            "          border-bottom-right-radius: 25px; }")
        report.add_css_rule(".fastqs {\n"
                            " border: 1px solid grey;\n"
                            " padding: 5px;\n"
                            " margin: 5px 20px;\n"
                            "}")
        report.add_css_rule(".fastq {\n"
                            " border: 2px solid lightgray;\n"
                            " padding: 5px;\n"
                            " margin: 5px;\n"
                            " float: left;\n"
                            "}")
        report.add_css_rule(".clear { clear: both; }")
        report.add_css_rule("table.metadata { margin: 10 10;\n"
                            "          border: solid 1px grey;\n"
                            "          background-color: white;\n"
                            "          font-size: 90%; }")
        report.add_css_rule("table.metadata tr td:first-child {\n"
                            "          background-color: grey;\n"
                            "          color: white;\n"
                            "          padding: 2px 5px;\n"
                            "          font-weight: bold; }")
        report.add_css_rule("table.summary { border: solid 1px grey;\n"
                            "                background-color: white;\n"
                            "                font-size: 80% }")
        report.add_css_rule("table.summary th { background-color: grey;\n"
                            "                   color: white;\n"
                            "                   padding: 2px 5px; }")
        report.add_css_rule("table.summary td { text-align: right; \n"
                            "                   padding: 2px 5px;\n"
                            "                   border-bottom: solid 1px lightgray; }")
        report.add_css_rule("table.fastq_summary tr td:first-child {\n"
                            "          background-color: grey;\n"
                            "          color: white;\n"
                            "          font-weight: bold; }")
        report.add_css_rule("table.fastq_summary tr td:first-child a {\n"
                            "          color: white;\n"
                            "          font-weight: bold; }")
        report.add_css_rule("table.fastqc_summary span.PASS { font-weight: bold;\n"
                            "                                 color: green; }")
        report.add_css_rule("table.fastqc_summary span.WARN { font-weight: bold;\n"
                            "                                 color: orange; }")
        report.add_css_rule("table.fastqc_summary span.FAIL { font-weight: bold;\n"
                            "                                 color: red; }")
        report.add_css_rule("table.programs th { text-align: left;\n"
                            "                    background-color: grey;\n"
                            "                    color: white;\n"
                            "                    padding: 2px 5px; }")
        report.add_css_rule("table.programs td { padding: 2px 5px;\n"
                            "border-bottom: solid 1px lightgray; }")
        report.add_css_rule("p { font-size: 85%;\n"
                            "    color: #808080; }")
        # Rules for printing
        report.add_css_rule("@media print\n"
                            "{\n"
                            "a { color: black; text-decoration: none; }\n"
                            ".sample { page-break-before: always; }\n"
                            "table th { border-bottom: solid 1px lightgray; }\n"
                            ".no_print { display: none; }\n"
                            "}")
        # Write the report
        report.write(filename)
        # Return the output filename
        return filename

class QCSample(object):
    """
    Class describing QC results for an AnalysisSample

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

    def verify(self,qc_dir):
        """
        Check QC products for this sample

        Checks that fastq_screens and FastQC files were found.
        Returns True if the QC products are present, False
        otherwise.
        """
        for fq_pair in self.fastq_pairs:
            if not fq_pair.verify(qc_dir):
                return False
        return True

class FastqSet(object):
    """
    Class describing a set of Fastq files

    A set can be a single or a pair of fastq files.

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
        Return list of Fastq files

        """
        return filter(lambda fq: fq is not None,
                      self._fastqs)

    def verify(self,qc_dir):
        """
        Check QC products for this Fastq pair

        Checks that fastq_screens and FastQC files were found.
        Returns True if the QC products are present, False
        otherwise.

        Arguments:
          qc_dir (str): path to the location of the QC
            output directory

        """
        illumina_qc = IlluminaQC()
        for fq in self._fastqs:
            if fq is None:
                continue
            present,missing = illumina_qc.check_outputs(fq,qc_dir)
            if missing:
                return False
        return True

class QCReport(Document):
    """
    """
    def __init__(self,project,title=None):
        """
        """
        # Store project
        self.project = project
        if title is None:
            title = "QC report: %s" % project.name
        # Initialise superclass
        Document.__init__(self,title)
        # Initialise tables
        self.metadata_table = self.init_metadata_table()
        self.summary_table = self.init_summary_table()
        # Initialise report sections
        self.preamble = self.init_preamble_section()
        self.summary = self.init_summary_section()
        # Add data
        self.report_metadata()
        for sample in self.project.samples:
            self.report_sample(sample)

    def init_metadata_table(self):
        """
        """
        metadata_tbl = Table(('item','value',))
        metadata_tbl.no_header()
        metadata_tbl.add_css_classes('metadata')
        return metadata_tbl

    def init_summary_table(self):
        """
        """
        if self.project.info.paired_end:
            fields = ('sample',
                      'fastqs',
                      'reads',
                      'fastqc_r1',
                      'boxplot_r1',
                      'screens_r1',
                      'fastqc_r2',
                      'boxplot_r2',
                      'screens_r2',)
        else:
            fields = ('sample',
                      'fastq',
                      'reads',
                      'fastqc_r1',
                      'boxplot_r1',
                      'screens_r1')
        field_descriptions = { 'sample': 'Sample',
                               'fastq' : 'Fastq',
                               'fastqs': 'Fastqs (R1/R2)',
                               'reads': '#reads',
                               'fastqc_r1': 'FastQC',
                               'boxplot_r1': 'Boxplot',
                               'screens_r1': 'Screens',
                               'fastqc_r2': 'FastQC',
                               'boxplot_r2': 'Boxplot',
                               'screens_r2': 'Screens' }
        summary_tbl = Table(fields,**field_descriptions)
        summary_tbl.add_css_classes('summary','fastq_summary')
        return summary_tbl

    def init_preamble_section(self):
        """
        """
        preamble = self.add_section()
        preamble.add("Report generated by auto_process %s on %s" %
                     (get_version(),time.asctime()))
        return preamble

    def init_summary_section(self):
        """
        """
        summary = self.add_section("Summary",name='summary')
        summary.add(self.metadata_table)
        summary.add("%d samples | %d fastqs" % (len(self.project.samples),
                                                len(self.project.fastqs)))
        summary.add(self.summary_table)
        return summary

    def report_metadata(self):
        """
        """
        metadata_items = ['user','PI','library_type','organism',]
        if self.project.info.single_cell_platform is not None:
            metadata_items.insert(3,'single_cell_platform')
        metadata_titles = {
            'user': 'User',
            'PI': 'PI',
            'library_type': 'Library type',
            'single_cell_platform': 'Single cell preparation platform',
            'organism': 'Organism',
        }
        for item in metadata_items:
            if self.project.info[item]:
                self.metadata_table.add_row(
                    item=metadata_titles[item],
                    value=self.project.info[item])

    def report_sample(self,sample):
        """
        """
        sample = QCSample(sample)
        # Create a new section for the sample
        sample_report = self.add_section(
            "Sample: %s" % sample.name,
            name="sample_%s" % sample.name)
        sample_report.add_css_classes('sample')
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
                                           self.project.qc_dir)
            fastq_pair.report(sample_report)
            # Add line in summary table
            if sample_name is not None:
                idx = self.summary_table.add_row(sample=Link(sample_name,
                                                             sample_report))
            else:
                idx = self.summary_table.add_row(sample="&nbsp;")
            fastq_pair.update_summary_table(self.summary_table,idx=idx)

class QCReportFastqPair(object):
    """
    """
    def __init__(self,fastqr1,fastqr2,qc_dir):
        """
        """
        self.fastqr1 = fastqr1
        self.fastqr2 = fastqr2
        self.qc_dir = qc_dir

    @property
    def paired_end(self):
        """
        """
        return (self.fastqr2 is not None)

    @property
    def r1(self):
        """
        """
        return QCReportFastq(self.fastqr1,self.qc_dir)

    @property
    def r2(self):
        """
        """
        if self.fastqr2 is not None:
            return QCReportFastq(self.fastqr2,self.qc_dir)
        return None

    def report(self,sample_report,attrs=None):
        """
        """
        # Attributes to report
        if attrs is None:
            attrs = ('fastqc','fastq_screen','program_versions')
        # Add container section for Fastq pair
        fastqs_report = sample_report.add_subsection()
        fastqs_report.add_css_classes('fastqs')
        # Create sections for individual Fastqs
        for fq in (self.r1,self.r2):
            if fq is None:
                continue
            fq_report = fastqs_report.add_subsection(fq.name,
                                                     name=fq.safe_name)
            fq_report.add_css_classes('fastq')
            # Add reports for each requested 'attribute'
            for attr in attrs:
                if attr == "fastqc":
                    # FastQC outputs
                    fq.report_fastqc(fq_report)
                elif attr == "fastq_screen":
                    # FastQScreens
                    fq.report_fastq_screens(fq_report)
                elif attr == "program_versions":
                    # Versions of programs used
                    new_section = fq.report_program_versions(fq_report)
                else:
                    raise KeyError("'%s': unrecognised reporting element "
                                   % attr)
        # Add an empty section to clear HTML floats
        clear = fastqs_report.add_subsection()
        clear.add_css_classes("clear")

    def update_summary_table(self,summary_table,idx=None,fields=None):
        """
        """
        # Fields to report
        if fields is None:
            if self.paired_end:
                fields = ('fastqs',
                          'reads',
                          'boxplot_r1','boxplot_r2',
                          'fastqc_r1','fastqc_r2',
                          'screens_r1','screens_r2')
            else:
                fields = ('fastq',
                          'reads',
                          'boxplot_r1',
                          'fastqc_r1',
                          'screens_r1')
        # Add row to summary table
        if idx is None:
            idx = summary_table.add_row()
        # Populate with data
        for field in fields:
            if field == "fastq":
                summary_table.set_value(idx,'fastq',
                                        Link(self.r1.name,
                                             "#%s" % self.r1.safe_name))
            elif field == "fastqs":
                summary_table.set_value(idx,'fastqs',
                                        "%s<br />%s" %
                                        (Link(self.r1.name,
                                              "#%s" % self.r1.safe_name),
                                         Link(self.r2.name,
                                              "#%s" % self.r2.safe_name)))
            elif field == "reads":
                summary_table.set_value(idx,'reads',
                                        pretty_print_reads(
                                            self.r1.fastqc.data.basic_statistics(
                                                'Total Sequences')))
            elif field == "boxplot_r1":
                summary_table.set_value(idx,'boxplot_r1',
                                        Img(self.r1.uboxplot(),
                                            href="#boxplot_%s" %
                                            self.r1.safe_name))
            elif field == "boxplot_r2":
                summary_table.set_value(idx,'boxplot_r2',
                                        Img(self.r2.uboxplot(),
                                            href="#boxplot_%s" %
                                            self.r2.safe_name))
            elif field == "fastqc_r1":
                summary_table.set_value(idx,'fastqc_r1',
                                        Img(self.r1.ufastqcplot(),
                                            href="#fastqc_%s" %
                                            self.r1.safe_name))
            elif field == "fastqc_r2":
                summary_table.set_value(idx,'fastqc_r2',
                                        Img(self.r2.ufastqcplot(),
                                            href="#fastqc_%s" %
                                            self.r2.safe_name))
            elif field == "screens_r1":
                summary_table.set_value(idx,'screens_r1',
                                        Img(self.r1.uscreenplot(),
                                            href="#fastq_screens_%s" %
                                            self.r1.safe_name))
            elif field == "screens_r2":
                summary_table.set_value(idx,'screens_r2',
                                        Img(self.r2.uscreenplot(),
                                            href="#fastq_screens_%s" %
                                            self.r2.safe_name))
            else:
                raise KeyError("'%s': unrecognised field for summary "
                               "table" % field)

class QCReportFastq(object):
    """
    Provides interface to QC outputs for Fastq file

    Atributes:
      name: basename of the Fastq
      path: path to the Fastq
      fastqc: Fastqc instance
      fastq_screen.names: list of FastQScreen names
      fastq_screen.NAME.description
      fastq_screen.NAME.png: associated PNG file
      fastq_screen.NAME.txt: associated TXT file
    """
    def __init__(self,fastq,qc_dir):
        """
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
            self.fastq_screen['names'].append(name)
            self.fastq_screen[name] = AttributeDictionary()
            self.fastq_screen[name]["description"] = name.replace('_',' ').title()
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
        """
        screens_report = document.add_subsection("Screens",
                                                 name="fastq_screens_%s" %
                                                 self.safe_name)
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

    def uboxplot(self,inline=True):
        """
        """
        return uboxplot(self.fastqc.data.path,inline=inline)

    def ufastqcplot(self,inline=True):
        """
        """
        return ufastqcplot(self.fastqc.summary.path,inline=inline)

    def uscreenplot(self,inline=True):
        """
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
