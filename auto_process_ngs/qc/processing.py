#!/usr/bin/env python
#
# Library for QC reporting of processing stages

import os
from bcftbx.TabFile import TabFile
from ..analysis import split_sample_name
from ..docwriter import Document
from ..docwriter import Table
from ..docwriter import List
from ..docwriter import Link
from ..docwriter import Para
from ..docwriter import Img
from ..docwriter import WarningIcon
from .. import css_rules
from .reporting import pretty_print_reads
from .plots import ustackedbar

#######################################################################
# Functions
#######################################################################

def get_absolute_file_path(p,base=None):
    """
    Get absolute path for supplied path

    Arguments:
      p (str): path
      base (str): optional, base directory to
        use if p is relative

    Returns:
      String: absolute path for p.
    """
    if not os.path.isabs(p):
        if base is not None:
            p = os.path.join(os.path.abspath(base),p)
        else:
            p = os.path.abspath(p)
    return p

def report_processing_qc(analysis_dir,html_file):
    """
    Generate HTML report for processing statistics

    Arguments:
      analysis_dir (AutoProcess): AutoProcess instance for
        the directory to report the processing from
      html_file (str): destination path and file name for
        HTML report
    """
    # Initialise the HTML report
    processing_qc = Document("Processing report for %s" %
                             os.path.basename(analysis_dir.analysis_dir))
    processing_qc.add_css_rule(css_rules.QC_REPORT_CSS_RULES)
    processing_qc.add_css_rule("table { font-size: 80%;\n"
                               "        font-family: sans-serif; }")
    processing_qc.add_css_rule("td { text-align: right; }")
    processing_qc.add_css_rule("p.warning { padding: 5px;\n"
                               "            border: solid 1px red;\n"
                               "            background-color: F5BCA9;\n"
                               "            color: red;\n"
                               "            font-weight: bold;\n"
                               "            border-radius: 10px;\n"
                               "            display: inline-block; }")
    processing_qc.add_css_rule(".warnings { padding: 2px;\n"
                               "            border: solid 3px red;\n"
                               "            background-color: F5BCA9;\n"
                               "            color: red;\n"
                               "            font-weight: bold;\n"
                               "            margin: 10px;\n"
                               "            border-radius: 10px;\n"
                               "            display: inline-block; }")
    processing_qc.add_css_rule("img { vertical-align: middle; }")
    processing_qc.add_css_rule(".hide { display: none; }")
    # Add table of contents
    toc = processing_qc.add_section("Contents",name="toc")
    toc_list = List()
    toc.add(toc_list)
    # Add warnings section
    # This will be hidden if there are no issues
    status = True
    warnings = processing_qc.add_section(css_classes=("warnings",))
    warnings.add(Para(WarningIcon(size=50),
                      "There are issues with one or more lanes or samples"))
    # Per-lane statistics
    per_lane_stats_file = analysis_dir.params.per_lane_stats_file
    if per_lane_stats_file is None:
        per_lane_stats_file = "per_lane_statistics.info"
    per_lane_stats_file = get_absolute_file_path(per_lane_stats_file,
                                                 base=analysis_dir.analysis_dir)
    if os.path.exists(per_lane_stats_file):
        per_lane_stats = processing_qc.add_section(
            "Per-lane statistics",
            name="per_lane_stats")
        stats = TabFile(per_lane_stats_file,
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
        toc_list.add_item(Link("Per-lane statistics",per_lane_stats))
    # Per lane by sample statistics
    per_lane_sample_stats_file = get_absolute_file_path(
        "per_lane_sample_stats.info",
        base=analysis_dir.analysis_dir)
    if os.path.exists(per_lane_sample_stats_file):
        per_lane_sample_stats = processing_qc.add_section(
            "Per-lane statistics by sample",
            name="per_lane_sample_stats")
        lane_toc_list = List()
        per_lane_sample_stats.add(lane_toc_list)
        # Store the data for each lane
        lane_data = list()
        with open(per_lane_sample_stats_file,'r') as stats:
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
                s.add(Para(WarningIcon(),
                           "No samples reported for this lane",
                           css_classes=('warning',)))
                has_warnings = True
            elif min([d['nreads'] for d in data['samples']]) == 0:
                # There are samples with no reads
                s.add(Para(WarningIcon(),
                           "One or more samples with no reads",
                           css_classes=('warning',)))
                has_warnings = True
            # Add link to lane for lane ToC
            link = Link("Lane %d" % lane,s)
            if not has_warnings:
                lane_toc_list.add_item(link)
            else:
                lane_toc_list.add_item(WarningIcon(),link)
                status = False
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
        toc_list.add_item(Link("Per-lane statistics by sample",
                               per_lane_sample_stats),
                          lane_toc_list)
    # Per fastq statistics
    stats_file = get_absolute_file_path("statistics_full.info",
                                        base=analysis_dir.analysis_dir)
    if not os.path.exists(stats_file):
        if analysis_dir.params.stats_file is not None:
            stats_file = analysis_dir.params.stats_file
        else:
            stats_file = "statistics.info"
    stats_file = get_absolute_file_path(stats_file,
                                        base=analysis_dir.analysis_dir)
    if os.path.exists(stats_file):
        per_file_stats = processing_qc.add_section(
            "Per-file statistics by project",
            name="per_file_stats")
        project_toc_list = List()
        per_file_stats.add(project_toc_list)
        stats = TabFile(stats_file,first_line_is_header=True)
        projects = sorted(list(set([d['Project'] for d in stats])))
        lanes = filter(lambda c: c.startswith('L'),stats.header())
        sample = None
        for project in projects:
            # Get subset of lines for this project
            subset = sorted(filter(lambda d: d['Project'] == project,stats),
                            key=lambda l: split_sample_name(l['Sample']))
            # Work out which lanes are included
            subset_lanes = filter(lambda l:
                                  reduce(lambda x,y: x or bool(y),
                                         [d[l] for d in subset],
                                         False),
                                  lanes)
            # Add a new section for this project
            s = per_file_stats.add_subsection(
                "%s" % project,
                name="per_file_stats_%s" % project
            )
            # Check for problems
            has_warnings = False
            for line in subset:
                nreads = filter(lambda n: n != '',
                                [line[l] for l in subset_lanes])
                if not nreads or min(nreads) == 0:
                    s.add(Para(WarningIcon(),"One or more Fastqs with zero "
                               "read counts in one or lanes",
                               css_classes=('warning',)))
                    has_warnings = True
                    break
            # Add link to project from ToC
            link = Link("%s" % project,s)
            if not has_warnings:
                project_toc_list.add_item(link)
            else:
                project_toc_list.add_item(WarningIcon(),link)
                status = False
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
                nreads = filter(lambda n: n != '',
                                [line[l] for l in subset_lanes])
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
        toc_list.add_item(Link("Per-file statistics by project",
                               per_file_stats),
                          project_toc_list)
    # Set the visibility of the warning header
    if status:
        warnings.add_css_classes("hide")
    # Add an non-visible section that the publisher can
    # read to determine if there were problems
    s = processing_qc.add_section(name="status",css_classes=("hide",))
    s.add("Status: %s" % ('OK' if status else 'WARNINGS',))
    # Write the processing QC summary file
    processing_qc.write(html_file)

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
