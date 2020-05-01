#!/usr/bin/env python
#
#     reporting: report processing QC from Fastq generation
#     Copyright (C) University of Manchester 2020 Peter Briggs
#

"""
Provides classes and functions for reporting QC information
about Fastq generation, based on data in files produced by the
``fastq_statistics.py`` utility.

The main class ``ProcessingQCReport`` takes these statistics
files as input and composes an HTML report.

The ``detect_processing_qc_warnings`` function scans an HTML
report from ``ProcessingQCReport`` and indicates whether the
report contains any warnings.
"""

#######################################################################
# Imports
#######################################################################

import os
import logging
from bcftbx.TabFile import TabFile
from ..analysis import split_sample_name
from ..docwriter import Document
from ..docwriter import Img
from ..docwriter import Link
from ..docwriter import List
from ..docwriter import Para
from ..docwriter import Table
from ..docwriter import WarningIcon
from ..qc.plots import ustackedbar
from ..qc.reporting import pretty_print_reads
from .. import css_rules

# Module specific logger
logger = logging.getLogger(__name__)

#######################################################################
# Classes
#######################################################################

class ProcessingQCReport(Document):
    """
    Create a processing QC report document

    Example usage:

    >>> report = ProcessingQCReport(...)
    >>> report.write("processing_qc_report.html")
    """
    def __init__(self,analysis_dir,stats_file,per_lane_stats_file,
                 per_lane_sample_stats_file):
        """
        Create a new ProcessingQCReport instance

        Arguments:
          analysis_dir (str): analysis directory to report
            the processing QC for
          stats_file (str): path to the 'statistics_full.info'
            statistics file
          per_lane_stats_file (str): path to the
            'per_lane_statistics.info' statistics file
          per_lane_sample_stats_file (str): path to the
            'per_lane_sample_stats.info' statistics file
        """
        # Initialise HTML report
        Document.__init__(self,"Processing report for %s" %
                          os.path.basename(analysis_dir))
        # Store locations of data files
        self._stats_file = stats_file
        self._per_lane_stats_file = per_lane_stats_file
        self._per_lane_sample_stats_file = per_lane_sample_stats_file
        # Add CSS styling rules
        self.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
        self.add_css_rule("table { font-size: 80%;\n"
                          "        font-family: sans-serif; }")
        self.add_css_rule("td { text-align: right; }")
        self.add_css_rule("p.warning { padding: 5px;\n"
                          "            border: solid 1px red;\n"
                          "            background-color: F5BCA9;\n"
                          "            color: red;\n"
                          "            font-weight: bold;\n"
                          "            border-radius: 10px;\n"
                          "            display: inline-block; }")
        self.add_css_rule(".warnings { padding: 2px;\n"
                          "            border: solid 3px red;\n"
                          "            background-color: F5BCA9;\n"
                          "            color: red;\n"
                          "            font-weight: bold;\n"
                          "            margin: 10px;\n"
                          "            border-radius: 10px;\n"
                          "            display: inline-block; }")
        self.add_css_rule("img { vertical-align: middle; }")
        self.add_css_rule(".hide { display: none; }")
        # Add warnings section
        self.init_warnings()
        # Add table of contents
        self.init_toc()
        # Add statistics sections
        self.add_per_lane_statistics()
        self.add_per_lane_sample_statistics()
        self.add_per_fastq_statistics()
        # Finalise the warnings
        self.finalise_warnings()

    def init_toc(self,title="Contents"):
        """
        Initialise the table of contents

        Arguments:
          title (str): title for the table of
            contents (default: "Contents")
        """
        toc = self.add_section(title,name="toc")
        self._toc_list = List()
        toc.add(self._toc_list)

    def add_to_toc(self,title,s,sub_toc=None):
        """
        Add an item to the table of contents

        Arguments:
          title (str): title for the item in the
            table of contents
          s (Object): document item to target in
            the table of contents
          sub_toc (List): optional, sub table of
            contents to append for this item
        """
        items = [Link(title,s)]
        if sub_toc is not None:
            items.append(sub_toc)
        self._toc_list.add_item(*items)

    def init_warnings(self):
        """
        Initialise section to indicate report has warnings

        The 'warnings' section is hidden by default, but a
        call to the 'flag_warnings' method will cause it
        to be unhidden when the HTML report is generated
        """
        self._warning_status = False
        self._warnings = self.add_section()
        self._warnings.add(
            self.warning("There are issues with one or more "
                         "lanes or samples",size=50))

    def flag_warnings(self):
        """
        Set flag to indicate that report includes warnings

        By default the flag is set to ``False`` (i.e. no
        warnings in the report); calling this method one
        or more time sets the flag to ``True`` (i.e. at
        least one warning in the report).

        The flag is checked by the ``finalise_warnings``
        method when the report is completed.
        """
        if not self._warning_status:
            self._warning_status = True

    def finalise_warnings(self):
        """
        Sets up warning indications on report

        This method should be called when the report
        content has been assembled; it performs the final
        actions for displaying the general warning
        indicator etc.
        """
        # Set the visibility of the warning header
        if not self._warning_status:
            self._warnings.add_css_classes("hide")
        # Add an non-visible section that the publisher can
        # read to determine if there were problems
        s = self.add_section(name="status",css_classes=("hide",))
        s.add("Status: %s" % ('OK' if not self._warning_status
                              else 'WARNINGS',))

    def warning(self,message,size=None):
        """
        Return object for displaying a warning

        Returns an object which will render to an HTML
        warning message (essentially a paragraph with
        a warning icon and the supplied message text).

        Arguments:
          message (str): warning message to display
          size (int): optional, set the size for the
            warning icon
        """
        return Para(WarningIcon(size=size),
                    message,
                    css_classes=('warning',))

    def add_per_lane_statistics(self):
        """
        Add a section with the per-lane statistics
        """
        # Per-lane statistics
        if not os.path.exists(self._per_lane_stats_file):
            logger.debug("No per-lane statistics file found")
            return
        per_lane_stats = self.add_section(
            "Per-lane statistics",
            name="per_lane_stats")
        stats = TabFile(self._per_lane_stats_file,
                        first_line_is_header=True)
        tbl = Table(columns=stats.header())
        tbl.append_columns("Assigned/unassigned")
        for line in stats:
            n = tbl.add_row()
            for c in stats.header():
                if c in ("Total reads","Assigned reads","Unassigned reads"):
                    value = pretty_print_reads(line[c])
                else:
                    value = line[c]
                tbl.set_value(n,c,value)
            tbl.set_value(n,"Assigned/unassigned",
                          Img(ustackedbar((line["Assigned reads"],
                                           line["Unassigned reads"]),
                                          length=100,height=15,
                                          colors=('red','white'),
                                              inline=True)))
        per_lane_stats.add(tbl)
        self.add_to_toc("Per-lane statistics",per_lane_stats)

    def add_per_lane_sample_statistics(self):
        """
        Add a section with the per-lane sample statistics
        """
        # Per lane by sample statistics
        if not os.path.exists(self._per_lane_sample_stats_file):
            logger.debug("No per-lane sample statistics file found")
            return
        per_lane_sample_stats = self.add_section(
            "Per-lane statistics by sample",
            name="per_lane_sample_stats")
        lane_toc_list = List()
        per_lane_sample_stats.add(lane_toc_list)
        # Store the data for each lane
        lane_data = list()
        # Read in the data from the stats file
        with open(self._per_lane_sample_stats_file,'r') as stats:
            for line in stats:
                if line.startswith("Lane "):
                    lane = int(line.split(' ')[-1])
                    lane_data.append({ 'lane': lane,
                                       'total_reads': None,
                                       'samples': [] })
                elif line.startswith("Total reads = "):
                    total_reads = int(line.split('=')[-1].strip())
                    lane_data[-1]['total_reads'] = total_reads
                elif line.startswith("- "):
                    pname = line.split()[1].split('/')[0]
                    sname = line.split()[1].split('/')[1]
                    nreads = int(line.split()[2])
                    percreads = line.split()[3]
                    lane_data[-1]['samples'].append(
                        { 'pname': pname,
                          'sname': sname,
                          'nreads': nreads,
                          'percreads': percreads })
        # Create a section and table for each lane
        for data in lane_data:
            lane = data['lane']
            s = per_lane_sample_stats.add_subsection(
                "Lane %d" % lane,
                name="per_lane_sample_stats_lane%d" % lane
            )
            # Check for problems
            has_warnings = False
            if not data['samples']:
                # No samples reported
                s.add(self.warning("No samples reported for this lane"))
                has_warnings = True
            elif min([d['nreads'] for d in data['samples']]) == 0:
                # There are samples with no reads
                s.add(self.warning("One or more samples with no reads"))
                has_warnings = True
            # Add link to lane for lane ToC
            link = Link("Lane %d" % lane,s)
            if not has_warnings:
                lane_toc_list.add_item(link)
            else:
                lane_toc_list.add_item(WarningIcon(),link)
                self.flag_warnings()
            # Write out the data, if there is any
            if not data['samples']:
                continue
            max_reads = max([d['nreads'] for d in data['samples']])
            total_reads = data['total_reads']
            current_project = None
            tbl = Table(columns=('pname','sname',
                                 'nreads','percreads',
                                 'barplot'),
                        pname='Project',
                        sname='Sample',
                        nreads='Nreads',
                        percreads='%reads',
                        barplot='',
                    )
            s.add(tbl)
            # Sort the sample data into order of sample name
            samples = sorted([s for s in data['samples']],
                             key=lambda s: split_sample_name(s['sname']))
            # Write the table
            for sample in samples:
                pname = sample['pname']
                sname = sample['sname']
                nreads = sample['nreads']
                percreads = sample['percreads']
                if pname == current_project:
                    pname = "&nbsp;"
                else:
                    current_project = pname
                barplot = ustackedbar((nreads,max_reads-nreads),
                                      length=100,height=5,
                                      colors=('black','lightgrey'),
                                      bbox=False,
                                      inline=True)
                if nreads == 0:
                    sname = Para(WarningIcon(),sname)
                tbl.add_row(pname=pname,
                            sname=sname,
                            nreads=pretty_print_reads(nreads),
                            percreads=percreads,
                            barplot=Img(barplot))
            tbl.add_row(pname="Total reads for lane %d" % lane,
                        nreads=pretty_print_reads(total_reads))
        # Add link to section from main ToC
        self.add_to_toc("Per-lane statistics by sample",
                        per_lane_sample_stats,
                        lane_toc_list)

    def add_per_fastq_statistics(self):
        """
        Add a section with the per-Fastq statistics
        """
        # Per fastq statistics
        if not os.path.exists(self._stats_file):
            logger.debug("No per-Fastq statistics file found")
            return
        per_file_stats = self.add_section(
            "Per-file statistics by project",
            name="per_file_stats")
        project_toc_list = List()
        per_file_stats.add(project_toc_list)
        stats = TabFile(self._stats_file,first_line_is_header=True)
        projects = sorted(list(set([d['Project'] for d in stats])))
        lanes = [c for c in stats.header() if c.startswith('L')]
        sample = None
        for project in projects:
            # Get subset of lines for this project
            subset = sorted([d for d in stats if d['Project'] == project],
                            key=lambda l: split_sample_name(l['Sample']))
            # Determine which lanes this project appears in
            subset_lanes = []
            for l in lanes:
                for d in subset:
                    if d[l]:
                        subset_lanes.append(l)
                        break
            # Add a new section for this project
            s = per_file_stats.add_subsection(
                "%s" % project,
                name="per_file_stats_%s" % project
            )
            # Check for problems
            has_warnings = False
            for line in subset:
                nreads = [line[l] for l in subset_lanes if line[l] != '']
                if not nreads or min(nreads) == 0:
                    s.add(self.warning("One or more Fastqs with zero read "
                                       "counts in one or more lanes"))
                    has_warnings = True
                    break
            # Add link to project from ToC
            link = Link("%s" % project,s)
            if not has_warnings:
                project_toc_list.add_item(link)
            else:
                project_toc_list.add_item(WarningIcon(),link)
                self.flag_warnings()
            # Build the data of data
            tbl = Table(columns=('Sample','Fastq','Size'))
            if subset_lanes:
                tbl.append_columns(*subset_lanes)
            tbl.append_columns('Barplot','Nreads')
            s.add(tbl)
            for line in subset:
                if sample == line['Sample']:
                    sname = "&nbsp;"
                else:
                    sample = line['Sample']
                    sname = sample
                data = { 'Sample': sname,
                         'Fastq': line['Fastq'],
                         'Size': line['Size'],
                         'Nreads': (pretty_print_reads(line['Nreads'])
                                    if line['Nreads'] != '' else '')
                }
                for l in subset_lanes:
                    data[l] = (pretty_print_reads(line[l])
                               if line[l] != '' else '')
                nreads = [line[l] for l in subset_lanes if line[l] != '']
                if not nreads:
                    nreads = [0,]
                if min(nreads) == 0:
                    # Add warning icon to Fastq with no reads in
                    # at least one lane
                    data['Fastq'] = Para(WarningIcon(),data['Fastq'])
                barplot = ustackedbar(nreads,
                                      length=100,height=10,
                                      colors=('grey','lightgrey'),
                                      bbox=True,
                                      inline=True)
                data['Barplot'] = Img(barplot)
                tbl.add_row(**data)
        # Add to table of contents
        self.add_to_toc("Per-file statistics by project",
                        per_file_stats,
                        project_toc_list)

#######################################################################
# Functions
#######################################################################

def detect_processing_qc_warnings(html_file):
    """
    Look for warning text in processing_qc.html file

    Arguments:
      html_file (str): path to HTML report file

    Returns:
      Boolean: True if warnings were found, False if not.
    """
    with open(html_file) as fp:
        for line in fp:
            if "Status: WARNINGS" in line:
                return True
                break
    return False
