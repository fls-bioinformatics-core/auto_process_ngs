#!/usr/bin/env python
#
# illumina_qc library

import sys
import os
import logging
import time
from bcftbx.IlluminaData import IlluminaFastq
from bcftbx.IlluminaData import cmp_sample_names
from bcftbx.TabFile import TabFile
from bcftbx.qc.report import strip_ngs_extensions
from bcftbx.utils import extract_prefix
from bcftbx.utils import extract_index
from ..docwriter import Document
from ..docwriter import Table
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import Target
from .fastqc import Fastqc
from .fastq_screen import Fastqscreen
from .plots import uscreenplot
from .plots import ufastqcplot
from .plots import uboxplot
from .plots import encode_png
from .. import get_version

FASTQ_SCREENS = ('model_organisms',
                 'other_organisms',
                 'rRNA',)

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class QCReporter:
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
        self._qc_dir = self._project.qc_dir
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

    def verify(self):
        """
        Check that the QC outputs are correct

        Returns True if the QC appears to have run successfully,
        False if not.

        """
        verified = True
        for sample in self._samples:
            if not sample.verify(self._qc_dir):
                verified = False
        return verified

    def report(self,title=None,filename=None,relative_links=False):
        """
        Report the QC for the project

        Arguments:
          title (str): optional, specify title for the report
            (defaults to '<PROJECT_NAME>: QC report')
          filename (str): optional, specify path and name for
            the output report file (defaults to
            '<PROJECT_NAME>.qc_report.html')
          relative_links (boolean): optional, if set to True
            then use relative paths for links in the report
            (default is to use absolute paths)

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
        report = Document(title=title)
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
        report.add_css_rule(".fastq_r1, .fastq_r2 {\n"
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
        # Preamble
        preamble = report.add_section()
        preamble.add("Report generated by auto_process %s on %s" %
                     (get_version(),time.asctime()))
        # Set up summary section & metadata
        summary = report.add_section("Summary",name='summary')
        metadata_tbl = Table(('item','value',))
        metadata_titles = {
            'user': 'User',
            'PI': 'PI',
            'library_type': 'Library type',
            'organism': 'Organism',
        }
        for item in ('user','PI','library_type','organism',):
            if self._project.info[item]:
                metadata_tbl.add_row(item=metadata_titles[item],
                                     value=self._project.info[item])
        if metadata_tbl.nrows > 0:
            metadata_tbl.no_header()
            metadata_tbl.add_css_classes('metadata')
            summary.add(metadata_tbl)
        # Add summary table for samples, fastqs and microplots
        summary.add("%d samples | %d fastqs" % (len(self._samples),
                                                len(self._project.fastqs)))
        summary_tbl = Table(('sample',),sample='Sample')
        summary_tbl.add_css_classes('summary','fastq_summary')
        summary.add(summary_tbl)
        if self.paired_end:
            summary_tbl.append_columns('fastqs',fastqs='Fastqs (R1/R2)')
        else:
            summary_tbl.append_columns('fastq',fastq='Fastq')
        summary_tbl.append_columns('reads','fastqc_r1','boxplot_r1','screens_r1',
                                   reads='#reads',fastqc_r1='FastQC',boxplot_r1='Boxplot',
                                   screens_r1='Screens')
        if self.paired_end:
            summary_tbl.append_columns('fastqc_r2','boxplot_r2','screens_r2',
                                   fastqc_r2='FastQC',boxplot_r2='Boxplot',
                                   screens_r2='Screens')
        # Write entries for samples, fastqs etc
        current_sample = None
        for i,sample in enumerate(self._samples):
            logger.debug("Reporting sample #%3d: %s " % (i+1,sample.name))
            sample_name = sample.name
            sample_report = report.add_section("Sample: %s" % sample_name,
                                               name="sample_%s" % sample_name)
            sample_report.add_css_classes('sample')
            if self.paired_end:
                sample_report.add("%d fastq R1/R2 pairs" %
                                  len(sample.fastq_pairs))
            else:
                sample_report.add("%d fastqs" % len(sample.fastq_pairs))
            for fq_pair in sample.fastq_pairs:
                # Sample name for first pair only
                if sample_name is not None:
                    idx = summary_tbl.add_row(sample=Link(sample_name,
                                                          sample_report))
                else:
                    idx = summary_tbl.add_row(sample="&nbsp;")
                # Container for fastqs
                fqs_report = sample_report.add_subsection()
                fqs_report.add_css_classes('fastqs')
                # Fastq name(s)
                fq_r1 = os.path.basename(fq_pair.r1)
                if self.paired_end:
                    fq_r2 = os.path.basename(fq_pair.r2)
                else:
                    fq_r2 = None
                if self.paired_end:
                    # Create subsections for R1 and R2
                    fqr1_report = fqs_report.add_subsection(fq_r1)
                    fqr2_report = fqs_report.add_subsection(fq_r2)
                    # Add classes
                    fqr1_report.add_css_classes('fastq_r1')
                    fqr2_report.add_css_classes('fastq_r2')
                    # Add entries to summary table
                    summary_tbl.set_value(idx,'fastqs',
                                          "%s<br />%s" %
                                          (Link(fq_r1,fqr1_report),
                                           Link(fq_r2,fqr2_report)))
                else:
                    # Create subsection for R1 only
                    fqr1_report = fqs_report.add_subsection(fq_r1)
                    # Add classes
                    fqr1_report.add_css_classes('fastq_r1')
                    # Add entry to summary table
                    summary_tbl.set_value(idx,'fastq',Link(fq_r1,
                                                           fqr1_report))
                self._report_fastq(fq_r1,'r1',summary_tbl,idx,
                                   fqr1_report,relpath=relpath)
                if self.paired_end:
                    self._report_fastq(fq_r2,'r2',summary_tbl,idx,
                                       fqr2_report,relpath=relpath)
                # Reset sample name for remaining pairs
                sample_name = None
                # Add an empty section to clear HTML floats
                clear = fqs_report.add_subsection()
                clear.add_css_classes("clear")
        # Write the report
        report.write(filename)

    def _report_fastq(self,fq,read_id,summary,idx,report,
                      relpath=None):
        """
        Generate report section for a Fastq file

        Arguments:
          fq (str): Fastq file that is being reported
          read_id (str): either 'r1' or 'r2'
          summary (Table): summary table object
          idx (integer): row index for summary table
          report (Section): container for the report

        """
        # Report FastQC results
        fastqc_report = report.add_subsection("FastQC")
        try:
            # Locate FastQC outputs
            fastqc = Fastqc(os.path.join(self._qc_dir,fastqc_output(fq)[0]))
            # FastQC quality boxplot
            fastqc_report.add("Per base sequence quality boxplot:")
            boxplot = Img(fastqc.quality_boxplot(inline=True),
                          height=250,
                          width=480,
                          href=fastqc.summary.link_to_module(
                              'Per base sequence quality',
                              relpath=relpath),
                          name="boxplot_%s" % fq)
            fastqc_report.add(boxplot)
            try:
                summary.set_value(idx,'boxplot_%s' % read_id,
                                  Img(uboxplot(fastqc.data.path,inline=True),
                                      href=boxplot))
            except Exception,ex:
                logger.error("Failed to generate boxplot for %s: %s"
                             % (fq,ex))
            # FastQC summary table
            fastqc_report.add("FastQC summary:")
            fastqc_tbl = Target("fastqc_%s" % fq)
            fastqc_report.add(fastqc_tbl,
                              fastqc.summary.html_table(relpath=relpath))
            if relpath:
                fastqc_html_report = os.path.relpath(fastqc.html_report,relpath)
            else:
                fastqc_html_report = fastqc.html_report
            fastqc_report.add("%s for %s" % (Link("Full FastQC report",
                                                  fastqc_html_report),
                                             fq))
            # Populate line in main Fastqs summary table
            if read_id == 'r1':
                nreads = fastqc.data.basic_statistics('Total Sequences')
                summary.set_value(idx,'reads',pretty_print_reads(nreads))
            try:
                summary.set_value(idx,'fastqc_%s' % read_id,
                                  Img(ufastqcplot(
                                      fastqc.summary.path,inline=True),
                                      href=fastqc_tbl))
            except Exception,ex:
                logger.error("Failed to generate Fastqc microplot for %s: %s"
                             % (fq,ex))
        except Exception,ex:
            # Unable to get the FastQC data
            logger.warning("Unable to load FastQC data for %s: %s" %
                           (fq,ex))
            # Add placeholders for missing data
            if read_id == 'r1':
                summary.set_value(idx,'reads','-')
            fastqc_report.add("!!!No FastQC data available!!!")
        # Report fastq_screen outputs
        screens_report = report.add_subsection("Screens")
        fastq_screens = Target("fastq_screens_%s" % fq)
        screens_report.add(fastq_screens)
        screen_files = []
        fastq_screen_txt = []
        for name in FASTQ_SCREENS:
            description = name.replace('_',' ').title()
            png,txt = fastq_screen_output(fq,name)
            png = os.path.join(self._qc_dir,png)
            txt = os.path.join(self._qc_dir,txt)
            if relpath:
                png_href = os.path.relpath(png,relpath)
                txt_href = os.path.relpath(txt,relpath)
            else:
                png_href = png
                txt_href = txt
            screens_report.add(description)
            if os.path.exists(png):
                screens_report.add(Img(encode_png(png),
                                       height=250,
                                       href=png_href))
            else:
                logger.warning("Unable to find screen PNG: %s" % png)
                screens_report.add("!!!No FastqScreen plot available!!!")
            if os.path.exists(txt):
                screen_files.append(txt)
                fastq_screen_txt.append(
                    Link(description,txt_href).html())
            else:
                logger.warning("Unable to find raw screen data: %s" % txt)
                screens_report.add("!!!No FastqScreen data available!!!")
        screens_report.add("Raw screen data: " +
                           " | ".join(fastq_screen_txt))
        try:
            summary.set_value(idx,'screens_%s' % read_id,
                              Img(uscreenplot(screen_files,inline=True),
                                  href=fastq_screens))
        except Exception,ex:
            logger.error("Failed to generate microscreen plots for %s: %s"
                         % (fq,ex))
        # Program versions
        versions = report.add_subsection("Program versions")
        versions.add(self._program_versions(fq))

    def _program_versions(self,fastq):
        """
        """
        # Program versions table
        tbl = Table(("Program","Version"))
        tbl.add_css_classes("programs","summary")
        # Fetch the version info
        try:
            fastqc_version = Fastqc(
                os.path.join(self._qc_dir,
                             fastqc_output(fastq)[0])).version
        except Exception,ex:
            logger.error("Unable to get Fastqc version for %s: %s"
                         % (fastq,ex))
            fastqc_version = "?"
        try:
            fastq_screen_version = Fastqscreen(
                os.path.join(self._qc_dir,
                             fastq_screen_output(fastq,
                                                 FASTQ_SCREENS[0])[1])).version
        except Exception,ex:
            logger.error("Unable to get Fastq_screen version for %s: %s"
                         % (fastq,ex))
            fastq_screen_version = "?"
        # Add to table
        tbl.add_row(Program='fastqc',Version=fastqc_version)
        tbl.add_row(Program='fastq_screen',Version=fastq_screen_version)
        return tbl

class QCSample:
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

class FastqSet:
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
        for fq in self._fastqs:
            if fq is None:
                continue
            present,missing = check_qc_outputs(fq,qc_dir)
            if missing:
                return False
        return True

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

def fastq_screen_output(fastq,screen_name):
    """
    Generate name of fastq_screen output files

    Given a Fastq file name and a screen name, the outputs from
    fastq_screen will look like:

    - {FASTQ}_{SCREEN_NAME}_screen.png
    - {FASTQ}_{SCREEN_NAME}_screen.txt

    Arguments:
       fastq (str): name of Fastq file
       screen_name (str): name of screen

    Returns:
       tuple: fastq_screen output names (without leading path)

    """
    base_name = "%s_%s_screen" % (strip_ngs_extensions(os.path.basename(fastq)),
                                  str(screen_name))
    
    return (base_name+'.png',base_name+'.txt')

def fastqc_output(fastq):
    """
    Generate name of FastQC outputs

    Given a Fastq file name, the outputs from FastQC will look
    like:

    - {FASTQ}_fastqc/
    - {FASTQ}_fastqc.html
    - {FASTQ}_fastqc.zip

    Arguments:
       fastq (str): name of Fastq file

    Returns:
       tuple: FastQC outputs (without leading paths)

    """
    base_name = "%s_fastqc" % strip_ngs_extensions(os.path.basename(fastq))
    return (base_name,base_name+'.html',base_name+'.zip')

def expected_qc_outputs(fastq,qc_dir):
    """
    Return list of expected QC products for a FASTQ file

    Arguments:
      fastq (str): name of FASTQ file
      qc_dir (str): path to QC directory

    Returns:
      List: list of paths to the expected associated QC products

    """
    expected = []
    # FastQC outputs
    expected.extend([os.path.join(qc_dir,f)
                     for f in fastqc_output(fastq)])
    # Fastq_screen outputs
    for name in FASTQ_SCREENS:
        expected.extend([os.path.join(qc_dir,f)
                         for f in fastq_screen_output(fastq,name)])
    return expected

def check_qc_outputs(fastq,qc_dir):
    """
    Return lists of present and missing QC products for FASTQ file

    Arguments:
      fastq (str): name of FASTQ file
      qc_dir (str): path to QC directory

    Returns:
      Tuple: tuple of the form (present,missing) where present,
        missing are lists of paths to associated QC products which
        are present in the QC dir, or are missing.

    """
    present = []
    missing = []
    # Check that outputs exist
    for output in expected_qc_outputs(fastq,qc_dir):
        if os.path.exists(output):
            present.append(output)
        else:
            missing.append(output)
    return (present,missing)

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
